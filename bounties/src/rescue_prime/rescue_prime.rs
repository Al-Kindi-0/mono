use ff::PrimeField;
use std::sync::Arc;

use super::rescue_prime_params::RescuePrimeParams;

#[derive(Clone, Debug)]
pub struct RescuePrime<S: PrimeField> {
    pub(crate) params: Arc<RescuePrimeParams<S>>,
}

impl<S: PrimeField> RescuePrime<S> {
    pub fn new(params: &Arc<RescuePrimeParams<S>>) -> Self {
        RescuePrime {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[S]) -> Vec<S> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();

        for r in 0..self.params.rounds {
            current_state = self.sbox(&current_state);
            current_state = self.affine(&current_state, 2 * r);
            current_state = self.sbox_inverse(&current_state);
            current_state = self.affine(&current_state, 2 * r + 1);
        }

        current_state
    }

    fn sbox(&self, input: &[S]) -> Vec<S> {
        input
            .iter()
            .map(|el| {
                let mut el2 = *el;
                el2.square();
                let res = match self.params.d {
                    3 => {
                        let mut out = el2;
                        out.mul_assign(&el);
                        out
                    }
                    5 => {
                        let mut out = el2;
                        out.square();
                        out.mul_assign(&el);
                        out
                    }
                    _ => {
                        assert!(false);
                        *el
                    }
                };
                res
            })
            .collect()
    }

    fn sbox_inverse(&self, input: &[S]) -> Vec<S> {
        input.iter().map(|el| el.pow(&self.params.d_inv)).collect()
    }

    fn affine(&self, input: &[S], round: usize) -> Vec<S> {
        let mat_result = self.matmul(input, &self.params.mds);
        Self::add_rc(&mat_result, &self.params.round_constants[round])
    }

    fn matmul(&self, input: &[S], mat: &[Vec<S>]) -> Vec<S> {
        let t = mat.len();
        debug_assert!(t == input.len());
        let mut out = vec![S::zero(); t];
        for row in 0..t {
            for col in 0..t {
                let mut tmp = mat[row][col];
                tmp.mul_assign(&input[col]);
                out[row].add_assign(&tmp);
            }
        }
        out
    }

    fn add_rc(input: &[S], round_constants: &[S]) -> Vec<S> {
        debug_assert!(input.len() == round_constants.len());
        input
            .iter()
            .zip(round_constants.iter())
            .map(|(a, b)| {
                let mut r = *a;
                r.add_assign(b);
                r
            })
            .collect()
    }
}

// #############################################################################

#[cfg(test)]
mod rescue_prime_kats {
    use super::*;

    use crate::fields::{field::Fp, utils};
    use crate::rescue_prime::rescue_prime_instances::*;

    use ff::{from_hex, Field};

    type Scalar = Fp;

    #[test]
    fn easy1_kats() {
        let rescue = RescuePrime::new(&RESCUE_PRIME_PARAMS_EASY1);
        let input: Vec<Scalar> = vec![Scalar::zero(), Scalar::one(), utils::from_u64::<Scalar>(2)];
        let perm = rescue.permutation(&input);
        assert_eq!(perm[0], from_hex("0x9775230f09921974").unwrap());
        assert_eq!(perm[1], from_hex("0x04b677d2493fe1d5").unwrap(),);
        assert_eq!(perm[2], from_hex("0x4015ef9a4c3afd18").unwrap(),);
    }

    #[test]
    fn easy2_kats() {
        let rescue = RescuePrime::new(&RESCUE_PRIME_PARAMS_EASY2);
        let input: Vec<Scalar> = vec![Scalar::zero(), Scalar::one(), utils::from_u64::<Scalar>(2)];
        let perm = rescue.permutation(&input);
        assert_eq!(perm[0], from_hex("0x469a02719e3ad1a4").unwrap());
        assert_eq!(perm[1], from_hex("0xd0e7e76477e414da").unwrap(),);
        assert_eq!(perm[2], from_hex("0xfcd1cb836abf8799").unwrap(),);
    }

    #[test]
    fn medium_kats() {
        let rescue = RescuePrime::new(&RESCUE_PRIME_PARAMS_MEDIUM);
        let input: Vec<Scalar> = vec![Scalar::zero(), Scalar::one(), utils::from_u64::<Scalar>(2)];
        let perm = rescue.permutation(&input);
        assert_eq!(perm[0], from_hex("0xa51fd5a1b81ee239").unwrap());
        assert_eq!(perm[1], from_hex("0x6f86f6ad1129cfc9").unwrap(),);
        assert_eq!(perm[2], from_hex("0xa93b46cb1310031e").unwrap(),);
    }

    #[test]
    fn hard1_kats() {
        let rescue = RescuePrime::new(&RESCUE_PRIME_PARAMS_HARD1);
        let input: Vec<Scalar> = vec![Scalar::zero(), Scalar::one(), utils::from_u64::<Scalar>(2)];
        let perm = rescue.permutation(&input);
        assert_eq!(perm[0], from_hex("0xb0aa37a1a00626a2").unwrap());
        assert_eq!(perm[1], from_hex("0xb8865d92a8bb736f").unwrap(),);
        assert_eq!(perm[2], from_hex("0xa5b5a8b00c10395a").unwrap(),);
    }

    #[test]
    fn hard2_kats() {
        let rescue = RescuePrime::new(&RESCUE_PRIME_PARAMS_HARD2);
        let input: Vec<Scalar> = vec![Scalar::zero(), Scalar::one(), utils::from_u64::<Scalar>(2)];
        let perm = rescue.permutation(&input);
        assert_eq!(perm[0], from_hex("0x03890b141a2ea2c6").unwrap());
        assert_eq!(perm[1], from_hex("0x562d54deec652cb3").unwrap(),);
        assert_eq!(perm[2], from_hex("0x8d7be2b0a0f3d4bf").unwrap(),);
    }
}
