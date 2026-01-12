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
    pub fn read_chunk(&self, chunk_idx: usize) -> std::io::Result<[stark_rings_linalg::SparseMatrix<R>; 3]> {
        use std::io::{BufReader, Seek, SeekFrom};
        const IO_BUFFER_SIZE: usize = 256 * 1024 * 1024;

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

        let mut chunk_matrices: [stark_rings_linalg::SparseMatrix<R>; 3] = std::array::from_fn(|_| {
            stark_rings_linalg::SparseMatrix {
                nrows,
                ncols: self.ncols,
                coeffs: Vec::with_capacity(nrows),
            }
        });

        for matrix in &mut chunk_matrices {
            for _ in 0..nrows {
                r.read_exact(&mut buf4)?;
                let num_terms = u32::from_le_bytes(buf4) as usize;
                let mut row = Vec::with_capacity(num_terms);
                for _ in 0..num_terms {
                    r.read_exact(&mut buf4)?;
                    let col_idx = u32::from_le_bytes(buf4) as usize;
                    r.read_exact(&mut buf8)?;
                    let coeff = i64::from_le_bytes(buf8);
                    if coeff == 0 {
                        continue;
                    }
                    let abs = coeff.unsigned_abs();
                    let base = <R as PolyRing>::BaseRing::from(abs);
                    let mut val = R::from(base);
                    if coeff < 0 {
                        val = -val;
                    }
                    row.push((val, col_idx));
                }
                matrix.coeffs.push(row);
            }
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
            return Ok(cache);
        }
    }

    // Build cache: read chunks from the .r1lf file and write them in a fast random-access format.
    // IMPORTANT: we must preserve signed i64 coefficients; do NOT roundtrip through `R`.
    let direct = R1LfChunkReader::open(path, chunk_size)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let num_chunks = direct.num_chunks();
    let ncols = expected_ncols;

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
    let mut src = std::fs::File::open(path)?;
    for i in 0..num_chunks {
        offsets[i] = w.stream_position()?;
        let start = i * chunk_size;
        let end = ((i + 1) * chunk_size).min(stats.num_constraints);
        let actual_rows = end - start;
        let padded_rows = next_power_of_two(actual_rows);

        // Chunk encoding: nrows (u64), then 3 matrices (A,B,C).
        w.write_all(&(padded_rows as u64).to_le_bytes())?;

        let (a0, b0, c0) = direct.chunk_offsets(i).map_err(|e| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, e)
        })?;
        write_matrix_chunk_from_r1lf(&mut src, &mut w, a0, actual_rows, padded_rows)?;
        write_matrix_chunk_from_r1lf(&mut src, &mut w, b0, actual_rows, padded_rows)?;
        write_matrix_chunk_from_r1lf(&mut src, &mut w, c0, actual_rows, padded_rows)?;
    }
    w.flush()?;

    // Backfill offsets.
    w.seek(SeekFrom::Start(offsets_pos))?;
    for off in offsets {
        w.write_all(&off.to_le_bytes())?;
    }
    w.flush()?;

    open_chunk_cache::<R>(&cache_path, &stats.digest)
}

fn open_chunk_cache<R>(path: &str, expected_digest: &[u8; 32]) -> std::io::Result<R1LfChunkCache<R>>
where
    R: OverField,
    R::BaseRing: Zq + From<u64>,
{
    use std::io::{BufReader, Read};
    const IO_BUFFER_SIZE: usize = 256 * 1024 * 1024;

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
    src: &mut std::fs::File,
    dst: &mut std::io::BufWriter<std::fs::File>,
    start_offset: u64,
    actual_rows: usize,
    padded_rows: usize,
) -> std::io::Result<()> {
    src.seek(SeekFrom::Start(start_offset))?;
    let mut buf4 = [0u8; 4];
    let mut buf12 = [0u8; 12];

    for _ in 0..actual_rows {
        src.read_exact(&mut buf4)?;
        let num_terms = u32::from_le_bytes(buf4) as usize;
        dst.write_all(&(num_terms as u32).to_le_bytes())?;
        for _ in 0..num_terms {
            src.read_exact(&mut buf12)?;
            // term is (u32 idx, i64 coeff)
            dst.write_all(&buf12)?;
        }
    }
    for _ in actual_rows..padded_rows {
        dst.write_all(&0u32.to_le_bytes())?;
    }
    Ok(())
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

fn scan_matrix_offsets(
    file: &mut File,
    start_offset: u64,
    num_constraints: usize,
    chunk_size: usize,
) -> Result<(Vec<u64>, u64), String> {
    file.seek(SeekFrom::Start(start_offset))
        .map_err(|e| format!("{e}"))?;

    let mut offsets = Vec::with_capacity((num_constraints + chunk_size - 1) / chunk_size);
    let mut pos = start_offset;

    for row_idx in 0..num_constraints {
        if row_idx % chunk_size == 0 {
            offsets.push(pos);
        }
        // num_terms: u32
        let mut buf4 = [0u8; 4];
        file.read_exact(&mut buf4).map_err(|e| format!("{e}"))?;
        pos += 4;
        let num_terms = u32::from_le_bytes(buf4) as usize;

        // skip terms: (u32 idx + i64 coeff) = 12 bytes each
        let skip = (num_terms as u64) * 12;
        file.seek(SeekFrom::Current(skip as i64))
            .map_err(|e| format!("{e}"))?;
        pos += skip;
    }

    Ok((offsets, pos))
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

