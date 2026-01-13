//! Loader for SP1 "R1LF" (LF-targeted lifted R1CS) files.
//!
//! This is a minimal, research-focused parser intended to bridge SP1 → LF+ without
//! going through Symphony's SP1 chunk-cache format.
//!
//! File format is written by `sp1_recursion_compiler::r1cs::lf::R1CSLf`.
//! It stores signed i64 coefficients (so we can represent `p_bb`).
//!
//! We support two modes:
//! - **Direct reader**: reads chunks by seeking into the `.r1lf` file.
//! - **Chunk cache**: Symphony-style `{path}.chunks` cache for fast random access and
//!   stable padded dimensions. This is the recommended path for LF+ experiments.

#![cfg(feature = "we_gate")]

use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};

use ark_ff::PrimeField;
use stark_rings::{OverField, PolyRing, Zq};

use crate::we_statement::WeParams;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Round up to next power of 2.
fn next_power_of_two(n: usize) -> usize {
    if n == 0 {
        return 1;
    }
    1usize << (usize::BITS - (n - 1).leading_zeros())
}

/// Compute the padded witness column count used by the R1LF chunk cache builder.
///
/// This is the canonical `ncols` to use when deriving WE parameters (e.g. `nvars_cm`),
/// and it intentionally includes any aux vars introduced by the lift (since `num_vars`
/// in the R1LF header already includes them).
pub fn padded_ncols_from_header(header: &R1LfHeader, pad_cols_to_multiple_of: usize) -> Result<usize, String> {
    if pad_cols_to_multiple_of == 0 {
        return Err("pad_cols_to_multiple_of must be > 0".to_string());
    }
    let blocks = (header.num_vars + pad_cols_to_multiple_of - 1) / pad_cols_to_multiple_of;
    let blocks_pow2 = next_power_of_two(blocks);
    Ok(blocks_pow2 * pad_cols_to_multiple_of)
}

/// Return \(nvars = log2(ncols)\) for power-of-two `ncols`.
pub fn nvars_from_ncols_pow2(ncols: usize) -> Result<usize, String> {
    if ncols == 0 || !ncols.is_power_of_two() {
        return Err(format!("ncols must be a power of two (got {ncols})"));
    }
    Ok(usize::BITS as usize - 1 - ncols.leading_zeros() as usize)
}

/// Default WE-gate parameters for SP1 BabyBear-in-Frog integration over an R1LF cache.
///
/// This is the **canonical** parameterization we want statement-bound in the WE arithmetization:
/// - boundedness digit base: `b = 2^16`
/// - digits per value: `k = 2`
///
/// `mlen` is the number of matrices being committed in the surrounding protocol layer.
pub fn sp1_default_we_params_for_r1lf_cache<R: PolyRing>(
    cache: &R1LfChunkCache<R>,
    kappa: u64,
    mlen: u64,
) -> Result<WeParams, String>
where
    R::BaseRing: PrimeField,
{
    let nvars = nvars_from_ncols_pow2(cache.ncols)?;
    // log_{d'}(q) where d' = d/2
    let lnq = (R::BaseRing::MODULUS_BIT_SIZE as f64) * std::f64::consts::LN_2;
    let l = (lnq / ((R::dimension() / 2) as f64).ln()).ceil() as u64;

    Ok(WeParams {
        nvars_setchk: nvars as u64,
        degree_setchk: 3,
        nvars_cm: nvars as u64,
        degree_cm: 2,
        kappa,
        ring_dim_d: R::dimension() as u64,
        decomp_b: 1u64 << 16,
        k: 2,
        l,
        mlen,
    })
}

/// Metadata parsed from the R1LF header.
#[derive(Debug, Clone)]
pub struct R1LfHeader {
    pub digest: [u8; 32],
    pub p_bb: u64,
    pub num_vars: usize,
    pub num_constraints: usize,
    pub num_public: usize,
    pub total_nonzeros: u64,
}

// ============================================================================
// Symphony-style chunk cache for R1LF
// ============================================================================

const R1LF_CHUNK_MAGIC: &[u8; 4] = b"LFC1"; // LF Chunk v1
const R1LF_CHUNK_VERSION: u32 = 1;
const R1LF_CHUNK_ZSTD_MAGIC: u32 = u32::from_le_bytes(*b"ZST1");
// Footer for O(1) cache integrity checks (seek-to-end).
// Old caches without this footer will be treated as incomplete and rebuilt.
const R1LF_CHUNK_FOOTER_MAGIC: u32 = u32::from_le_bytes(*b"LFF1");
const R1LF_CHUNK_FOOTER_LEN: u64 = 4 + 32 + 8; // magic + digest + file_len

fn write_chunk_cache_footer(
    w: &mut dyn std::io::Write,
    digest: &[u8; 32],
    file_len: u64,
) -> std::io::Result<()> {
    w.write_all(&R1LF_CHUNK_FOOTER_MAGIC.to_le_bytes())?;
    w.write_all(digest)?;
    w.write_all(&file_len.to_le_bytes())?;
    Ok(())
}

/// Random-access reader for a `{path}.chunks` cache file (loads one chunk at a time).
pub struct R1LfChunkCache<R> {
    pub stats: R1LfHeader,
    pub chunk_size: usize,
    pub ncols: usize,
    pub num_chunks: usize,
    pub(crate) cache_path: String,
    pub(crate) chunk_offsets: Vec<u64>, // absolute file offsets
    _phantom: std::marker::PhantomData<R>,
}

impl<R> R1LfChunkCache<R>
where
    R: OverField + PolyRing,
    R::BaseRing: Zq + PrimeField + From<u64>,
{
    /// Read one chunk as **base-ring** sparse matrices.
    ///
    /// SP1 R1LF coefficients are const-coeff in the cyclotomic ring model, so materializing full
    /// `R` elements per nonzero is wasted work and memory. This returns matrices over `R::BaseRing`
    /// directly, which is what downstream SP1/LF+ checks and streaming code paths use.
    pub fn read_chunk(
        &self,
        chunk_idx: usize,
    ) -> std::io::Result<[stark_rings_linalg::SparseMatrix<<R as PolyRing>::BaseRing>; 3]> {
        use std::io::{BufReader, Seek, SeekFrom};
        // This reader sits on top of a file slice for a *single* matrix payload; keep it modest.
        // (Huge buffers here can dominate runtime via page-zeroing.)
        const IO_BUFFER_SIZE: usize = 8 * 1024 * 1024;

        if chunk_idx >= self.num_chunks {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "chunk_idx out of range",
            ));
        }
        let file = std::fs::File::open(&self.cache_path)?;
        let mut r = BufReader::with_capacity(IO_BUFFER_SIZE, file);
        r.seek(SeekFrom::Start(self.chunk_offsets[chunk_idx]))?;

        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];
        r.read_exact(&mut buf8)?;
        let nrows = u64::from_le_bytes(buf8) as usize;
        r.read_exact(&mut buf4)?;
        let magic = u32::from_le_bytes(buf4);
        if magic != R1LF_CHUNK_ZSTD_MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "R1LF chunk cache: unsupported chunk payload (expected ZST1)",
            ));
        }

        // Actual rows for this chunk (last chunk may be partial).
        let start = chunk_idx * self.chunk_size;
        let end = ((chunk_idx + 1) * self.chunk_size).min(self.stats.num_constraints);
        let actual_rows = end.saturating_sub(start);

        let mut chunk_matrices: [stark_rings_linalg::SparseMatrix<<R as PolyRing>::BaseRing>; 3] =
            std::array::from_fn(|_| {
            stark_rings_linalg::SparseMatrix {
                nrows,
                ncols: self.ncols,
                coeffs: Vec::with_capacity(nrows),
            }
        });

        for matrix in &mut chunk_matrices {
            r.read_exact(&mut buf8)?;
            let clen = u64::from_le_bytes(buf8) as u64;
            let mut limited = r.by_ref().take(clen);
            let dec = zstd::stream::Decoder::new(&mut limited)?;
            let mut dec = std::io::BufReader::with_capacity(64 * 1024 * 1024, dec);

            // Decode only the actual rows, then pad with empty rows to `nrows`.
            for _ in 0..actual_rows {
                dec.read_exact(&mut buf4)?;
                let num_terms = u32::from_le_bytes(buf4) as usize;
                let mut row = Vec::with_capacity(num_terms);
                for _ in 0..num_terms {
                    dec.read_exact(&mut buf4)?;
                    let col_idx = u32::from_le_bytes(buf4) as usize;
                    dec.read_exact(&mut buf8)?;
                    let coeff = i64::from_le_bytes(buf8);
                    if coeff == 0 {
                        continue;
                    }
                    let abs = coeff.unsigned_abs();
                    let mut val = <R as PolyRing>::BaseRing::from(abs);
                    if coeff < 0 {
                        val = -val;
                    }
                    row.push((val, col_idx));
                }
                matrix.coeffs.push(row);
            }
            for _ in actual_rows..nrows {
                matrix.coeffs.push(Vec::new());
            }

            // Ensure we consume the full compressed payload so the underlying reader advances.
            std::io::copy(&mut dec, &mut std::io::sink())?;
            // Drain any remaining compressed bytes from `limited` (should be empty).
            std::io::copy(&mut limited, &mut std::io::sink())?;
        }
        Ok(chunk_matrices)
    }
}

/// Open (or build) the `{path}.chunks` cache and return a random-access reader.
pub fn open_sp1_r1lf_chunk_cache<R>(
    path: &str,
    chunk_size: usize,
    pad_cols_to_multiple_of: usize,
) -> std::io::Result<R1LfChunkCache<R>>
where
    R: OverField + PolyRing + Clone + Send + Sync,
    R::BaseRing: Zq + PrimeField + From<u64> + Send + Sync,
{
    #[derive(Clone, Copy, Debug)]
    struct ChunkMeta {
        a0: u64,
        a1: u64,
        b0: u64,
        b1: u64,
        c0: u64,
        c1: u64,
        actual_rows: usize,
        padded_rows: usize,
    }

    if pad_cols_to_multiple_of == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "pad_cols_to_multiple_of must be > 0",
        ));
    }

    let stats = read_r1lf_stats(path)?;
    let cache_path = format!("{path}.chunks");

    let blocks = (stats.num_vars + pad_cols_to_multiple_of - 1) / pad_cols_to_multiple_of;
    let blocks_pow2 = next_power_of_two(blocks);
    let expected_ncols = blocks_pow2 * pad_cols_to_multiple_of;

    if let Ok(cache) = open_chunk_cache::<R>(&cache_path, &stats.digest) {
        if cache.chunk_size == chunk_size && cache.ncols == expected_ncols {
            // Cache exists; sanity-check that it isn't truncated/corrupt.
            if cache_file_seems_complete(&cache)? {
                return Ok(cache);
            }
            // Otherwise rebuild.
            let _ = std::fs::remove_file(&cache_path);
        }
    }

    // Build cache: read chunks from the .r1lf file and write them in a fast random-access format.
    // IMPORTANT: we must preserve signed i64 coefficients; do NOT roundtrip through `R`.
    let direct = R1LfChunkReader::open(path, chunk_size)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let num_chunks = direct.num_chunks();
    let ncols = expected_ncols;

    // Precompute source offsets per chunk (sequential; cheap) so parallel workers don't need to share `direct`.
    let mut metas: Vec<ChunkMeta> = Vec::with_capacity(num_chunks);
    let (_a_start, b_start, c_start) = direct
        .chunk_offsets(0)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let file_len = std::fs::metadata(path)?.len();
    for i in 0..num_chunks {
        let start = i * chunk_size;
        let end = ((i + 1) * chunk_size).min(stats.num_constraints);
        let actual_rows = end - start;
        let padded_rows = next_power_of_two(actual_rows);
        let (a0, b0, c0) = direct
            .chunk_offsets(i)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let (a1, b1, c1) = if i + 1 < num_chunks {
            direct
                .chunk_offsets(i + 1)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?
        } else {
            (b_start, c_start, file_len)
        };
        metas.push(ChunkMeta {
            a0,
            a1,
            b0,
            b1,
            c0,
            c1,
            actual_rows,
            padded_rows,
        });
    }

    // Parallel cache build (feature = "parallel"): build per-chunk blobs in parallel, then concatenate.
    #[cfg(feature = "parallel")]
    {
        use std::io::BufWriter;

        let tmp_dir = format!("{cache_path}.tmpd");
        let _ = std::fs::remove_dir_all(&tmp_dir);
        std::fs::create_dir_all(&tmp_dir)?;

        let src_path = path.to_string();

        metas
            .par_iter()
            .enumerate()
            .try_for_each(|(i, meta)| -> std::io::Result<()> {
                // Keep per-thread buffers small to avoid blowing RAM with many Rayon threads.
                const SRC_BUF: usize = 8 * 1024 * 1024;
                const DST_BUF: usize = 8 * 1024 * 1024;

                let src_file = std::fs::File::open(&src_path)?;
                let mut src = std::io::BufReader::with_capacity(SRC_BUF, src_file);

                let chunk_path = format!("{tmp_dir}/{i}.chunk");
                let dst_file = std::fs::File::create(&chunk_path)?;
                let mut dst = BufWriter::with_capacity(DST_BUF, dst_file);

                // Chunk encoding:
                //   - nrows (u64) = padded_rows
                //   - magic (u32) = ZST1
                //   - for each of 3 matrices: (clen:u64, zstd(bytes))
                dst.write_all(&(meta.padded_rows as u64).to_le_bytes())?;
                dst.write_all(&R1LF_CHUNK_ZSTD_MAGIC.to_le_bytes())?;
                // Must flush before any `write_matrix_chunk_from_r1lf` backfills lengths via `seek`.
                dst.flush()?;
                write_matrix_chunk_from_r1lf(&mut src, &mut dst, meta.a0, meta.a1, meta.actual_rows, meta.padded_rows)?;
                write_matrix_chunk_from_r1lf(&mut src, &mut dst, meta.b0, meta.b1, meta.actual_rows, meta.padded_rows)?;
                write_matrix_chunk_from_r1lf(&mut src, &mut dst, meta.c0, meta.c1, meta.actual_rows, meta.padded_rows)?;
                dst.flush()?;
                Ok(())
            })?;

        // Concatenate deterministically into the final cache file.
        let file = std::fs::File::create(&cache_path)?;
        let mut w = std::io::BufWriter::with_capacity(256 * 1024 * 1024, file);

        // Header (fixed)
        w.write_all(R1LF_CHUNK_MAGIC)?;
        w.write_all(&R1LF_CHUNK_VERSION.to_le_bytes())?;
        w.write_all(&stats.digest)?;
        w.write_all(&stats.p_bb.to_le_bytes())?;
        w.write_all(&(stats.num_vars as u64).to_le_bytes())?;
        w.write_all(&(stats.num_constraints as u64).to_le_bytes())?;
        w.write_all(&(stats.num_public as u64).to_le_bytes())?;
        w.write_all(&stats.total_nonzeros.to_le_bytes())?;
        w.write_all(&(chunk_size as u64).to_le_bytes())?;
        w.write_all(&(ncols as u64).to_le_bytes())?;
        w.write_all(&(num_chunks as u64).to_le_bytes())?;

        // Offset table: backfilled after writing chunks.
        let offsets_pos = w.stream_position()?;
        for _ in 0..num_chunks {
            w.write_all(&0u64.to_le_bytes())?;
        }
        w.flush()?;

        let mut offsets = vec![0u64; num_chunks];
        for i in 0..num_chunks {
            offsets[i] = w.stream_position()?;
            let chunk_path = format!("{tmp_dir}/{i}.chunk");
            let mut r = std::io::BufReader::with_capacity(
                64 * 1024 * 1024,
                std::fs::File::open(&chunk_path)?,
            );
            std::io::copy(&mut r, &mut w)?;
        }
        w.flush()?;

        // Backfill offsets.
        w.seek(SeekFrom::Start(offsets_pos))?;
        for off in offsets {
            w.write_all(&off.to_le_bytes())?;
        }
        w.flush()?;

        // Append footer: O(1) integrity check on open (seek-to-end).
        let end_pos = w.seek(SeekFrom::End(0))?;
        write_chunk_cache_footer(&mut w, &stats.digest, end_pos + R1LF_CHUNK_FOOTER_LEN)?;
        w.flush()?;

        // Best-effort cleanup of temp chunks.
        let _ = std::fs::remove_dir_all(&tmp_dir);
    }

    // Single-threaded cache build (no `parallel` feature).
    #[cfg(not(feature = "parallel"))]
    {
        let file = std::fs::File::create(&cache_path)?;
        let mut w = std::io::BufWriter::with_capacity(256 * 1024 * 1024, file);

        // Header (fixed)
        w.write_all(R1LF_CHUNK_MAGIC)?;
        w.write_all(&R1LF_CHUNK_VERSION.to_le_bytes())?;
        w.write_all(&stats.digest)?;
        w.write_all(&stats.p_bb.to_le_bytes())?;
        w.write_all(&(stats.num_vars as u64).to_le_bytes())?;
        w.write_all(&(stats.num_constraints as u64).to_le_bytes())?;
        w.write_all(&(stats.num_public as u64).to_le_bytes())?;
        w.write_all(&stats.total_nonzeros.to_le_bytes())?;
        w.write_all(&(chunk_size as u64).to_le_bytes())?;
        w.write_all(&(ncols as u64).to_le_bytes())?;
        w.write_all(&(num_chunks as u64).to_le_bytes())?;

        // Offset table: backfilled after writing chunks.
        let offsets_pos = w.stream_position()?;
        for _ in 0..num_chunks {
            w.write_all(&0u64.to_le_bytes())?;
        }
        w.flush()?;

        let mut offsets = vec![0u64; num_chunks];
        // Use a large buffered reader for the source `.r1lf` to avoid tiny read syscalls
        // (`read_exact` of 4/12 bytes) dominating CPU time during cache build.
        let src_file = std::fs::File::open(path)?;
        let mut src = std::io::BufReader::with_capacity(256 * 1024 * 1024, src_file);
        for (i, meta) in metas.iter().enumerate() {
            offsets[i] = w.stream_position()?;
            // Chunk encoding:
            //   - nrows (u64) = padded_rows
            //   - magic (u32) = ZST1
            //   - for each of 3 matrices: (clen:u64, zstd(bytes))
            w.write_all(&(meta.padded_rows as u64).to_le_bytes())?;
            w.write_all(&R1LF_CHUNK_ZSTD_MAGIC.to_le_bytes())?;
            // Must flush before any `write_matrix_chunk_from_r1lf` backfills lengths via `seek`.
            w.flush()?;
            write_matrix_chunk_from_r1lf(&mut src, &mut w, meta.a0, meta.a1, meta.actual_rows, meta.padded_rows)?;
            write_matrix_chunk_from_r1lf(&mut src, &mut w, meta.b0, meta.b1, meta.actual_rows, meta.padded_rows)?;
            write_matrix_chunk_from_r1lf(&mut src, &mut w, meta.c0, meta.c1, meta.actual_rows, meta.padded_rows)?;
        }
        w.flush()?;

        // Backfill offsets.
        w.seek(SeekFrom::Start(offsets_pos))?;
        for off in offsets {
            w.write_all(&off.to_le_bytes())?;
        }
        w.flush()?;

        // Append footer: O(1) integrity check on open (seek-to-end).
        let end_pos = w.seek(SeekFrom::End(0))?;
        write_chunk_cache_footer(&mut w, &stats.digest, end_pos + R1LF_CHUNK_FOOTER_LEN)?;
        w.flush()?;
    }

    open_chunk_cache::<R>(&cache_path, &stats.digest)
}

fn cache_file_seems_complete<R>(cache: &R1LfChunkCache<R>) -> std::io::Result<bool> {
    // Quick structural checks that catch the common failure mode: cache file exists and has a valid
    // header/digest, but the chunk payload is truncated (e.g. killed mid-build).
    use std::io::{Read, Seek, SeekFrom};

    let meta = std::fs::metadata(&cache.cache_path)?;
    let file_len = meta.len();
    // Must have at least header + footer.
    if file_len < (4 + 4 + 32 + 8 * 6 + 8 * 3) as u64 + R1LF_CHUNK_FOOTER_LEN {
        return Ok(false);
    }
    if cache.chunk_offsets.is_empty() {
        return Ok(false);
    }

    // Offsets should be monotone and within file.
    let mut prev = 0u64;
    for &off in &cache.chunk_offsets {
        if off < prev || off >= file_len {
            return Ok(false);
        }
        prev = off;
    }

    // O(1) integrity check: verify a fixed-size footer at end-of-file.
    let mut f = std::fs::File::open(&cache.cache_path)?;
    f.seek(SeekFrom::End(-(R1LF_CHUNK_FOOTER_LEN as i64)))?;
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];
    let mut digest = [0u8; 32];
    f.read_exact(&mut buf4)?;
    let magic = u32::from_le_bytes(buf4);
    if magic != R1LF_CHUNK_FOOTER_MAGIC {
        return Ok(false);
    }
    f.read_exact(&mut digest)?;
    if digest != cache.stats.digest {
        return Ok(false);
    }
    f.read_exact(&mut buf8)?;
    let claimed_len = u64::from_le_bytes(buf8);
    if claimed_len != file_len {
        return Ok(false);
    }
    Ok(true)
}

fn open_chunk_cache<R>(path: &str, expected_digest: &[u8; 32]) -> std::io::Result<R1LfChunkCache<R>>
where
    R: OverField,
    R::BaseRing: Zq + From<u64>,
{
    use std::io::{BufReader, Read};
    // We only read a small fixed header + an offsets table here; large buffers just add latency.
    const IO_BUFFER_SIZE: usize = 256 * 1024;

    let file = std::fs::File::open(path)?;
    let mut r = BufReader::with_capacity(IO_BUFFER_SIZE, file);

    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != R1LF_CHUNK_MAGIC {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Invalid R1LF chunk cache magic",
        ));
    }
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];
    r.read_exact(&mut buf4)?;
    let version = u32::from_le_bytes(buf4);
    if version != R1LF_CHUNK_VERSION {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "R1LF chunk cache version mismatch",
        ));
    }
    let mut digest = [0u8; 32];
    r.read_exact(&mut digest)?;
    if &digest != expected_digest {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "R1LF chunk cache digest mismatch",
        ));
    }
    r.read_exact(&mut buf8)?;
    let p_bb = u64::from_le_bytes(buf8);
    r.read_exact(&mut buf8)?;
    let num_vars = u64::from_le_bytes(buf8) as usize;
    r.read_exact(&mut buf8)?;
    let num_constraints = u64::from_le_bytes(buf8) as usize;
    r.read_exact(&mut buf8)?;
    let num_public = u64::from_le_bytes(buf8) as usize;
    r.read_exact(&mut buf8)?;
    let total_nonzeros = u64::from_le_bytes(buf8);
    r.read_exact(&mut buf8)?;
    let chunk_size = u64::from_le_bytes(buf8) as usize;
    r.read_exact(&mut buf8)?;
    let ncols = u64::from_le_bytes(buf8) as usize;
    r.read_exact(&mut buf8)?;
    let num_chunks = u64::from_le_bytes(buf8) as usize;

    let mut offsets = vec![0u64; num_chunks];
    for i in 0..num_chunks {
        r.read_exact(&mut buf8)?;
        offsets[i] = u64::from_le_bytes(buf8);
    }

    Ok(R1LfChunkCache {
        stats: R1LfHeader { digest, p_bb, num_vars, num_constraints, num_public, total_nonzeros },
        chunk_size,
        ncols,
        num_chunks,
        cache_path: path.to_string(),
        chunk_offsets: offsets,
        _phantom: std::marker::PhantomData,
    })
}

/// Chunked reader for R1LF.
pub struct R1LfChunkReader {
    file: File,
    header: R1LfHeader,
    chunk_size: usize,
    // Byte offsets (from file start) to each chunk start, per matrix.
    a_offsets: Vec<u64>,
    b_offsets: Vec<u64>,
    c_offsets: Vec<u64>,
}

impl R1LfChunkReader {
    pub fn open(path: &str, chunk_size: usize) -> Result<Self, String> {
        let mut file = File::open(path).map_err(|e| format!("{e}"))?;
        let header = read_header(&mut file)?;

        // Try to load cached offsets (symphony-style "chunk cache", but lightweight).
        let idx_path = format!("{path}.idx");
        let (a_offsets, b_offsets, c_offsets) =
            match try_load_idx(&idx_path, &header, chunk_size) {
                Ok(Some((a, b, c))) => (a, b, c),
                Ok(None) => {
                    // Scan once to compute chunk offsets for A, then B, then C.
                    let (a, a_end) =
                        scan_matrix_offsets(&mut file, 80, header.num_constraints, chunk_size)?;
                    let (b, b_end) =
                        scan_matrix_offsets(&mut file, a_end, header.num_constraints, chunk_size)?;
                    let (c, _c_end) =
                        scan_matrix_offsets(&mut file, b_end, header.num_constraints, chunk_size)?;
                    // Best-effort write.
                    let _ = write_idx(&idx_path, &header, chunk_size, &a, &b, &c);
                    (a, b, c)
                }
                Err(_e) => {
                    // If idx is corrupt, rescan and overwrite.
                    let (a, a_end) =
                        scan_matrix_offsets(&mut file, 80, header.num_constraints, chunk_size)?;
                    let (b, b_end) =
                        scan_matrix_offsets(&mut file, a_end, header.num_constraints, chunk_size)?;
                    let (c, _c_end) =
                        scan_matrix_offsets(&mut file, b_end, header.num_constraints, chunk_size)?;
                    let _ = write_idx(&idx_path, &header, chunk_size, &a, &b, &c);
                    (a, b, c)
                }
            };

        Ok(Self { file, header, chunk_size, a_offsets, b_offsets, c_offsets })
    }

    #[inline]
    pub fn header(&self) -> &R1LfHeader {
        &self.header
    }

    #[inline]
    pub fn num_chunks(&self) -> usize {
        (self.header.num_constraints + self.chunk_size - 1) / self.chunk_size
    }

    pub fn read_chunk<R>(&mut self, chunk_idx: usize) -> Result<[stark_rings_linalg::SparseMatrix<R>; 3], String>
    where
        R: OverField + PolyRing,
        R::BaseRing: Zq + PrimeField + From<u64> + Send + Sync,
    {
        let num_chunks = self.num_chunks();
        if chunk_idx >= num_chunks {
            return Err(format!("chunk_idx out of range: {chunk_idx} (num_chunks={num_chunks})"));
        }

        let start_row = chunk_idx * self.chunk_size;
        let end_row = ((chunk_idx + 1) * self.chunk_size).min(self.header.num_constraints);
        let nrows = end_row - start_row;
        let ncols = self.header.num_vars;

        let a0 = *self
            .a_offsets
            .get(chunk_idx)
            .ok_or_else(|| "missing A offset".to_string())?;
        let b0 = *self
            .b_offsets
            .get(chunk_idx)
            .ok_or_else(|| "missing B offset".to_string())?;
        let c0 = *self
            .c_offsets
            .get(chunk_idx)
            .ok_or_else(|| "missing C offset".to_string())?;

        let a = read_matrix_chunk::<R>(&mut self.file, a0, nrows, ncols)?;
        let b = read_matrix_chunk::<R>(&mut self.file, b0, nrows, ncols)?;
        let c = read_matrix_chunk::<R>(&mut self.file, c0, nrows, ncols)?;

        Ok([a, b, c])
    }

    fn chunk_offsets(&self, chunk_idx: usize) -> Result<(u64, u64, u64), String> {
        let num_chunks = self.num_chunks();
        if chunk_idx >= num_chunks {
            return Err(format!("chunk_idx out of range: {chunk_idx} (num_chunks={num_chunks})"));
        }
        Ok((
            *self.a_offsets.get(chunk_idx).ok_or("missing A offset")?,
            *self.b_offsets.get(chunk_idx).ok_or("missing B offset")?,
            *self.c_offsets.get(chunk_idx).ok_or("missing C offset")?,
        ))
    }
}

fn read_header(file: &mut File) -> Result<R1LfHeader, String> {
    let mut hdr = [0u8; 80];
    file.read_exact(&mut hdr).map_err(|e| format!("{e}"))?;
    if &hdr[0..4] != b"R1LF" {
        return Err("Invalid R1LF magic".to_string());
    }
    let version = u32::from_le_bytes(hdr[4..8].try_into().unwrap());
    if version != 1 {
        return Err(format!("Unsupported R1LF version: {version}"));
    }
    let mut digest = [0u8; 32];
    digest.copy_from_slice(&hdr[8..40]);
    let p_bb = u64::from_le_bytes(hdr[40..48].try_into().unwrap());
    let num_vars = u64::from_le_bytes(hdr[48..56].try_into().unwrap()) as usize;
    let num_constraints = u64::from_le_bytes(hdr[56..64].try_into().unwrap()) as usize;
    let num_public = u64::from_le_bytes(hdr[64..72].try_into().unwrap()) as usize;
    let total_nonzeros = u64::from_le_bytes(hdr[72..80].try_into().unwrap());
    Ok(R1LfHeader { digest, p_bb, num_vars, num_constraints, num_public, total_nonzeros })
}

pub fn read_r1lf_stats(path: &str) -> std::io::Result<R1LfHeader> {
    let mut file = std::fs::File::open(path)?;
    read_header(&mut file).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

fn write_matrix_chunk_from_r1lf(
    src: &mut std::io::BufReader<std::fs::File>,
    dst: &mut std::io::BufWriter<std::fs::File>,
    start_offset: u64,
    end_offset: u64,
    actual_rows: usize,
    padded_rows: usize,
) -> std::io::Result<()> {
    // Compress the raw row encoding bytes for this chunk directly from the `.r1lf`.
    //
    // We keep the on-disk term encoding identical (per row: u32 num_terms, then `num_terms` terms
    // of 12 bytes), but store it compressed (zstd) to reduce disk footprint and I/O.
    if end_offset < start_offset {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "end_offset < start_offset",
        ));
    }
    dst.flush()?;
    let len = end_offset - start_offset;
    src.seek(SeekFrom::Start(start_offset))?;
    let mut limited = src.by_ref().take(len);

    // Placeholder for compressed length.
    let len_pos = dst.get_mut().stream_position()?;
    dst.write_all(&0u64.to_le_bytes())?;
    dst.flush()?;

    // Compress into `dst` and compute compressed length by file positions.
    let data_pos = dst.get_mut().stream_position()?;
    {
        // level=3: good speed/ratio tradeoff (default).
        let mut enc = zstd::stream::Encoder::new(dst.by_ref(), 3)?;
        std::io::copy(&mut limited, &mut enc)?;
        enc.finish()?;
    }
    dst.flush()?;
    let end_pos = dst.get_mut().stream_position()?;
    let clen = end_pos - data_pos;

    // Backfill compressed length.
    dst.flush()?;
    dst.get_mut().seek(SeekFrom::Start(len_pos))?;
    dst.get_mut().write_all(&(clen as u64).to_le_bytes())?;
    dst.get_mut().seek(SeekFrom::Start(end_pos))?;

    // We do NOT write padded rows into the compressed payload. Readers add empty rows up to
    // `padded_rows` (which is stored as the chunk nrows).
    let _ = (actual_rows, padded_rows);
    Ok(())
}

fn scan_matrix_offsets(
    file: &mut File,
    start_offset: u64,
    num_constraints: usize,
    chunk_size: usize,
) -> Result<(Vec<u64>, u64), String> {
    // Fast scanner for chunk start offsets.
    use std::io::{BufReader, Read, Seek, SeekFrom};
    const IO_BUFFER_SIZE: usize = 256 * 1024 * 1024;
    const DISCARD_BUF_SIZE: usize = 4 * 1024 * 1024;

    let mut f = file.try_clone().map_err(|e| format!("{e}"))?;
    f.seek(SeekFrom::Start(start_offset))
        .map_err(|e| format!("{e}"))?;
    let mut r = BufReader::with_capacity(IO_BUFFER_SIZE, f);

    let mut offsets = Vec::with_capacity((num_constraints + chunk_size - 1) / chunk_size);
    let mut pos = start_offset;
    let mut buf4 = [0u8; 4];
    let mut discard = vec![0u8; DISCARD_BUF_SIZE];

    #[inline]
    fn read_and_discard(
        r: &mut BufReader<std::fs::File>,
        mut n: usize,
        scratch: &mut [u8],
    ) -> Result<(), String> {
        use std::io::Read;
        while n > 0 {
            let take = n.min(scratch.len());
            r.read_exact(&mut scratch[..take])
                .map_err(|e| format!("{e}"))?;
            n -= take;
        }
        Ok(())
    }

    for row_idx in 0..num_constraints {
        if row_idx % chunk_size == 0 {
            offsets.push(pos);
        }
        r.read_exact(&mut buf4).map_err(|e| format!("{e}"))?;
        pos += 4;
        let num_terms = u32::from_le_bytes(buf4) as usize;
        let skip = num_terms * 12;
        read_and_discard(&mut r, skip, &mut discard)?;
        pos += skip as u64;
    }

    Ok((offsets, pos))
}

fn try_load_idx(
    idx_path: &str,
    header: &R1LfHeader,
    chunk_size: usize,
) -> Result<Option<(Vec<u64>, Vec<u64>, Vec<u64>)>, String> {
    let mut f = match File::open(idx_path) {
        Ok(f) => f,
        Err(_) => return Ok(None),
    };
    let mut hdr = [0u8; 56];
    f.read_exact(&mut hdr).map_err(|e| format!("{e}"))?;
    if &hdr[0..4] != b"R1LI" {
        return Ok(None);
    }
    let version = u32::from_le_bytes(hdr[4..8].try_into().unwrap());
    if version != 1 {
        return Ok(None);
    }
    let mut digest = [0u8; 32];
    digest.copy_from_slice(&hdr[8..40]);
    if digest != header.digest {
        return Ok(None);
    }
    let cs = u64::from_le_bytes(hdr[40..48].try_into().unwrap()) as usize;
    if cs != chunk_size {
        return Ok(None);
    }
    let num_chunks = u64::from_le_bytes(hdr[48..56].try_into().unwrap()) as usize;

    let mut read_u64_vec = |n: usize| -> Result<Vec<u64>, String> {
        let mut out = vec![0u64; n];
        let mut buf8 = [0u8; 8];
        for i in 0..n {
            f.read_exact(&mut buf8).map_err(|e| format!("{e}"))?;
            out[i] = u64::from_le_bytes(buf8);
        }
        Ok(out)
    };

    let a = read_u64_vec(num_chunks)?;
    let b = read_u64_vec(num_chunks)?;
    let c = read_u64_vec(num_chunks)?;
    Ok(Some((a, b, c)))
}

fn write_idx(
    idx_path: &str,
    header: &R1LfHeader,
    chunk_size: usize,
    a: &[u64],
    b: &[u64],
    c: &[u64],
) -> Result<(), String> {
    let mut f = File::create(idx_path).map_err(|e| format!("{e}"))?;
    f.write_all(b"R1LI").map_err(|e| format!("{e}"))?;
    f.write_all(&1u32.to_le_bytes()).map_err(|e| format!("{e}"))?;
    f.write_all(&header.digest).map_err(|e| format!("{e}"))?;
    f.write_all(&(chunk_size as u64).to_le_bytes())
        .map_err(|e| format!("{e}"))?;
    f.write_all(&(a.len() as u64).to_le_bytes())
        .map_err(|e| format!("{e}"))?;
    for vec in [a, b, c] {
        for &x in vec {
            f.write_all(&x.to_le_bytes()).map_err(|e| format!("{e}"))?;
        }
    }
    Ok(())
}

fn read_matrix_chunk<R>(
    file: &mut File,
    start_offset: u64,
    nrows: usize,
    ncols: usize,
) -> Result<stark_rings_linalg::SparseMatrix<R>, String>
where
    R: OverField + PolyRing,
    R::BaseRing: Zq + PrimeField + From<u64> + Send + Sync,
{
    file.seek(SeekFrom::Start(start_offset))
        .map_err(|e| format!("{e}"))?;

    let mut coeffs: Vec<Vec<(R, usize)>> = Vec::with_capacity(nrows);
    for _ in 0..nrows {
        let mut buf4 = [0u8; 4];
        file.read_exact(&mut buf4).map_err(|e| format!("{e}"))?;
        let num_terms = u32::from_le_bytes(buf4) as usize;
        let mut row: Vec<(R, usize)> = Vec::with_capacity(num_terms);

        for _ in 0..num_terms {
            let mut buf12 = [0u8; 12];
            file.read_exact(&mut buf12).map_err(|e| format!("{e}"))?;
            let idx = u32::from_le_bytes(buf12[0..4].try_into().unwrap()) as usize;
            let coeff = i64::from_le_bytes(buf12[4..12].try_into().unwrap());
            if coeff == 0 {
                continue;
            }
            let abs = coeff.unsigned_abs();
            let base = <R as PolyRing>::BaseRing::from(abs);
            let mut r = R::from(base);
            if coeff < 0 {
                r = -r;
            }
            row.push((r, idx));
        }
        coeffs.push(row);
    }

    Ok(stark_rings_linalg::SparseMatrix { nrows, ncols, coeffs })
}
