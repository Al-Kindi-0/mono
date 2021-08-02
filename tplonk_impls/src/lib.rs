//! # tplonk
//!
//! A pure Rust implementation of the ReinforcedConcrete Permutation
#![cfg_attr(feature = "asm", feature(asm))]

pub extern crate ff;

pub mod feistel_mimc;
pub mod fields;
pub mod pedersen_hash;
pub mod poseidon;
pub mod reinforced_concrete;
pub mod reinforced_concrete_st;
pub mod rescue;
pub mod rescue_prime;
pub mod sinsemilla;
