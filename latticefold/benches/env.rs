use std::{env, str::FromStr};

use lazy_static::lazy_static;

fn get_env_var<T: FromStr>(key: &str) -> Option<T> {
    if let Ok(var) = env::var(key) {
        if let Ok(v) = var.parse::<T>() {
            Some(v)
        } else {
            None
        }
    } else {
        None
    }
}

#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct Env {
    pub duration: f32,
    pub warmup: f32,
    pub GoldilocksRingNTT: bool,
    pub StarkRingNTT: bool,
    pub BabyBearRingNTT: bool,
    pub FrogRingNTT: bool,
    pub prover: bool,
    pub verifier: bool,
    pub ajtai: bool,
    pub linearization: bool,
    pub decomposition: bool,
    pub folding: bool,
    pub e2e: bool,
    pub x_len: Option<usize>,
    pub kappa: Option<usize>,
    pub w: Option<usize>,
    pub wit_len: Option<usize>,
    pub b: Option<u128>,
    pub l: Option<usize>,
    pub b_small: Option<usize>,
    pub k: Option<usize>,
}

lazy_static! {
    pub static ref ENV: Env = {
        let goldilocks = get_env_var::<String>("GOLDILOCKS").is_some();
        let stark = get_env_var::<String>("STARK").is_some();
        let babybear = get_env_var::<String>("BABYBEAR").is_some();
        let frog = get_env_var::<String>("FROG").is_some();

        let (goldilocks, stark, babybear, frog) = (
            goldilocks || !stark && !babybear && !frog,
            stark || !goldilocks && !babybear && !frog,
            babybear || !goldilocks && !stark && !frog,
            frog || !goldilocks && !stark && !babybear,
        );

        let prover = get_env_var::<String>("PROVER").is_some();
        let verifier = get_env_var::<String>("VERIFIER").is_some();
        let ajtai = get_env_var::<String>("AJTAI").is_some();
        let (prover, verifier, ajtai) = (
            prover || !verifier && !ajtai,
            verifier || !prover && !ajtai,
            ajtai || !prover && !verifier,
        );

        let linearization = get_env_var::<String>("LINEARIZATION").is_some();
        let decomposition = get_env_var::<String>("DECOMPOSITION").is_some();
        let folding = get_env_var::<String>("FOLDING").is_some();
        let e2e = get_env_var::<String>("E2E").is_some();

        let (linearization, decomposition, folding, e2e) = (
            linearization || !decomposition && !folding && !e2e,
            decomposition || !linearization && !folding && !e2e,
            folding || !linearization && !decomposition && !e2e,
            e2e || !linearization && !decomposition && !folding,
        );

        Env {
            duration: get_env_var("DURATION").unwrap_or(50.0),
            warmup: get_env_var("WARMUP").unwrap_or(1.0),
            GoldilocksRingNTT: goldilocks,
            StarkRingNTT: stark,
            BabyBearRingNTT: babybear,
            FrogRingNTT: frog,
            prover,
            verifier,
            ajtai,
            linearization,
            decomposition,
            folding,
            e2e,
            x_len: get_env_var("X_LEN"),
            kappa: get_env_var("KAPPA"),
            w: get_env_var("W"),
            wit_len: get_env_var("WIT_LEN"),
            b: get_env_var("B"),
            l: get_env_var("L"),
            b_small: get_env_var("B_SMALL"),
            k: get_env_var("K"),
        }
    };
}
