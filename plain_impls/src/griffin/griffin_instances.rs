use lazy_static::lazy_static;
use std::sync::Arc;

use crate::{
    fields::{bls12::FpBLS12, bn256::FpBN256, st::FpST},
    griffin::griffin_params::GriffinParams,
};

lazy_static! {
    // BN256
    pub static ref GRIFFIN_BN_PARAMS: Arc<GriffinParams<FpBN256>> = Arc::new(GriffinParams::new(3, 5, 12));
    // BLS12
    pub static ref GRIFFIN_BLS_PARAMS: Arc<GriffinParams<FpBLS12>> = Arc::new(GriffinParams::new(3, 5, 12));
    // ST
    pub static ref GRIFFIN_ST_PARAMS: Arc<GriffinParams<FpST>> = Arc::new(GriffinParams::new(3, 3, 16));
}
