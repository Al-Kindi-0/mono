use ff::{LegendreSymbol, PrimeField, SqrtField};

use crate::fields::utils;
use sha3::{digest::ExtendableOutput, digest::Update, Sha3XofReader, Shake128};

#[derive(Clone, Debug)]
pub struct ReinforcedConcreteParams<F: PrimeField + SqrtField> {
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) betas: [F; 2],
    pub(crate) si: Vec<u16>,
    pub(crate) divisor_i: Vec<u64>,
    pub(crate) reciprokal_i: Vec<u64>,
    pub(crate) norm_shift_i: Vec<u32>,
    pub(crate) sbox: Vec<u16>,
    pub(crate) d: usize,
}

impl<F: PrimeField + SqrtField> ReinforcedConcreteParams<F> {
    pub const PRE_ROUNDS: usize = 3;
    pub const POST_ROUNDS: usize = 3;
    pub const TOTAL_ROUNDS: usize = Self::PRE_ROUNDS + Self::POST_ROUNDS + 1;
    pub const T: usize = 3;
    pub const INIT_SHAKE: &'static str = "ReinforcedConcrete";

    pub fn new(d: usize, si: &[u16], sbox: &[u16]) -> Self {
        assert!(sbox.len() <= u16::MAX as usize);

        let mut shake = Self::init_shake();
        let betas = Self::instantiate_betas(&mut shake);
        let round_constants = Self::instantiate_rc(&mut shake);

        let len = si.len();
        let mut divisor_i = Vec::with_capacity(len);
        let mut reciprokal_i = Vec::with_capacity(len);
        let mut norm_shift_i = Vec::with_capacity(len);
        for s in si {
            let (div, rec) = utils::compute_normalized_divisor_and_reciproical(*s);
            divisor_i.push(div);
            reciprokal_i.push(rec);
            norm_shift_i.push((*s as u64).leading_zeros());
        }

        ReinforcedConcreteParams {
            round_constants,
            betas,
            si: si.to_owned(),
            divisor_i,
            reciprokal_i,
            norm_shift_i,
            sbox: Self::pad_sbox(sbox, si),
            d,
        }
    }

    fn init_shake() -> Sha3XofReader {
        let mut shake = Shake128::default();
        shake.update(Self::INIT_SHAKE);
        for i in F::char().as_ref() {
            shake.update(u64::to_le_bytes(*i));
        }
        shake.finalize_xof()
    }

    fn instantiate_betas(shake: &mut Sha3XofReader) -> [F; 2] {
        let mut betas = [F::zero(); 2];
        let one = F::one();
        let two = utils::from_u64::<F>(2);

        loop {
            let beta = utils::field_element_from_shake(shake);
            if beta != F::zero() && Self::is_quadratic_non_residue(&one, &beta) {
                betas[0] = beta;
                break;
            }
        }

        loop {
            let beta = utils::field_element_from_shake(shake);
            if beta != F::zero() && Self::is_quadratic_non_residue(&two, &beta) {
                betas[1] = beta;
                break;
            }
        }

        betas
    }

    fn is_quadratic_non_residue(alpha: &F, beta: &F) -> bool {
        let mut symbol = alpha.to_owned();
        symbol.square();
        let mut tmp = beta.to_owned();
        tmp.double();
        tmp.double();
        symbol.sub_assign(&tmp);
        symbol.legendre() == LegendreSymbol::QuadraticNonResidue
    }

    fn pad_sbox(sbox: &[u16], si: &[u16]) -> Vec<u16> {
        let len = sbox.len();

        let max = si.iter().max().expect("si are empty...").to_owned();
        let mut out = sbox.to_owned();

        out.reserve((max as usize) - len);
        for i in (len as u16)..max {
            out.push(i);
        }

        out
    }

    fn instantiate_rc(shake: &mut Sha3XofReader) -> Vec<Vec<F>> {
        (0..=Self::TOTAL_ROUNDS)
            .map(|_| {
                (0..Self::T)
                    .map(|_| utils::field_element_from_shake(shake))
                    .collect()
            })
            .collect()
    }

    pub fn get_t(&self) -> usize {
        Self::T
    }

    pub fn get_rounds(&self) -> usize {
        Self::TOTAL_ROUNDS
    }
}
