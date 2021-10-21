use lazy_static::lazy_static;
use std::sync::Arc;

use crate::{feistel_mimc::feistel_mimc_params::FeistelMimcParams, fields::field::Fp};

lazy_static! {
    pub static ref FM_PARAMS_EASY1: Arc<FeistelMimcParams<Fp>> =
        Arc::new(FeistelMimcParams::new(3, 6));
    pub static ref FM_PARAMS_EASY2: Arc<FeistelMimcParams<Fp>> =
        Arc::new(FeistelMimcParams::new(3, 10));
    pub static ref FM_PARAMS_MEDIUM: Arc<FeistelMimcParams<Fp>> =
        Arc::new(FeistelMimcParams::new(3, 14));
    pub static ref FM_PARAMS_HARD1: Arc<FeistelMimcParams<Fp>> =
        Arc::new(FeistelMimcParams::new(3, 18));
    pub static ref FM_PARAMS_HARD2: Arc<FeistelMimcParams<Fp>> =
        Arc::new(FeistelMimcParams::new(3, 22));
}