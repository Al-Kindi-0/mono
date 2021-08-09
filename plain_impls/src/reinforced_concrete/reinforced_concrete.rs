use std::sync::Arc;

use ff::{PrimeField, SqrtField};

use crate::fields::utils;

use super::reinforced_concrete_params::ReinforcedConcreteParams;

#[derive(Clone, Debug)]
pub struct ReinforcedConcrete<F: PrimeField + SqrtField> {
    pub(crate) params: Arc<ReinforcedConcreteParams<F>>,
}

impl<F: PrimeField + SqrtField> ReinforcedConcrete<F> {
    pub fn new(params: &Arc<ReinforcedConcreteParams<F>>) -> Self {
        debug_assert!(ReinforcedConcreteParams::<F>::T == 3);
        ReinforcedConcrete {
            params: Arc::clone(params),
        }
    }

    pub fn concrete(&self, state: &mut [F; 3], round: usize) {
        // multiplication by circ(2 1 1) is equal to state + sum(state)

        let mut sum = state[0];
        state.iter().skip(1).for_each(|el| sum.add_assign(el));

        for (el, rc) in state
            .iter_mut()
            .zip(self.params.round_constants[round].iter())
        {
            el.add_assign(&sum);
            el.add_assign(rc); // add round constant
        }
    }

    pub fn bricks(&self, state: &[F; 3]) -> [F; 3] {
        let mut new_state: [F; 3] = [F::zero(); 3];

        // squaring
        let mut x1_sq = state[0];
        x1_sq.square();
        let mut x2_sq = state[1];
        x2_sq.square();

        // x1
        let mut x1 = x1_sq;
        match self.params.d {
            3 => {}
            5 => x1.square(),
            _ => panic!("not implemented!"),
        }
        x1.mul_assign(&state[0]);
        new_state[0] = x1;

        // x2
        x1_sq.add_assign(&state[0]);
        x1_sq.add_assign(&self.params.betas[0]);
        x1_sq.mul_assign(&state[1]);
        new_state[1] = x1_sq;

        // x3
        x2_sq.add_assign(&state[1]);
        x2_sq.add_assign(&state[1]);
        x2_sq.add_assign(&self.params.betas[1]);
        x2_sq.mul_assign(&state[2]);
        new_state[2] = x2_sq;

        new_state
    }

    pub fn decompose(&self, val: &F) -> Vec<u16> {
        let len = self.params.si.len();
        let mut res = vec![0; len];
        let mut repr = val.into_repr();

        for i in (1..self.params.si.len()).rev() {
            let (r, m) = utils::divide_long_using_recip::<F>(
                &repr,
                self.params.divisor_i[i],
                self.params.reciprokal_i[i],
                self.params.norm_shift_i[i],
            );
            repr = r;
            res[i] = m;
        }

        res[0] = repr.as_ref()[0] as u16;

        // just debugging
        if cfg!(debug_assertions) {
            let repr_ref = repr.as_ref();
            debug_assert!(repr_ref[0] < self.params.si[0] as u64);
            repr_ref
                .iter()
                .skip(1)
                .for_each(|el| debug_assert!(*el == 0));
        }

        res
    }

    pub fn compose(&self, vals: &[u16]) -> F {
        let mut repr = F::Repr::default();
        repr.as_mut()[0] = vals[0] as u64;

        for (val, s) in vals.iter().zip(self.params.si.iter()).skip(1) {
            repr = utils::mul_by_single_word::<F>(&repr, *s as u64);
            repr = utils::add_single_word::<F>(&repr, *val as u64);
        }
        F::from_repr(repr).unwrap()
    }

    pub fn bars(&self, state: &[F; 3]) -> [F; 3] {
        let mut s = state.to_owned();
        for el in s.iter_mut() {
            let mut vals = self.decompose(&el);
            for val in vals.iter_mut() {
                // *val = self.params.sbox[*val as usize];
                // safe because sbox is padded to the correct size in params
                unsafe {
                    *val = *self.params.sbox.get_unchecked(*val as usize);
                }
            }
            *el = self.compose(&vals);
        }
        s
    }

    pub fn permutation(&self, input: &[F; 3]) -> [F; 3] {
        assert_eq!(ReinforcedConcreteParams::<F>::T, input.len());
        let mut current_state = input.to_owned();
        // first concrete
        self.concrete(&mut current_state, 0);

        // first rounds
        for i in 1..=ReinforcedConcreteParams::<F>::PRE_ROUNDS {
            current_state = self.bricks(&current_state);
            self.concrete(&mut current_state, i);
        }

        // bar round
        current_state = self.bars(&current_state);
        self.concrete(
            &mut current_state,
            ReinforcedConcreteParams::<F>::PRE_ROUNDS + 1,
        );

        // final rounds
        for i in ReinforcedConcreteParams::<F>::PRE_ROUNDS + 2
            ..=ReinforcedConcreteParams::<F>::TOTAL_ROUNDS
        {
            current_state = self.bricks(&current_state);
            self.concrete(&mut current_state, i);
        }
        current_state
    }

    pub fn hash(&self, el1: &F, el2: &F) -> F {
        let input: [F; 3] = [el1.to_owned(), el2.to_owned(), F::zero()];
        self.permutation(&input)[0]
    }
}

#[cfg(test)]
mod reinforced_concrete_tests_bn256 {
    use ff::{from_hex, Field};

    use crate::{
        fields::bn256::FpBN256, reinforced_concrete::reinforced_concrete_instances::RC_BN_PARAMS,
    };

    type Scalar = FpBN256;

    use super::*;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let rc = ReinforcedConcrete::new(&RC_BN_PARAMS);
        for _ in 0..TESTRUNS {
            let input1: [Scalar; 3] = [
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
            ];

            let mut input2: [Scalar; 3];
            loop {
                input2 = [
                    utils::random_scalar(true),
                    utils::random_scalar(true),
                    utils::random_scalar(true),
                ];
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = rc.permutation(&input1);
            let perm2 = rc.permutation(&input1);
            let perm3 = rc.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn compose() {
        let rc = ReinforcedConcrete::new(&RC_BN_PARAMS);

        for _ in 0..TESTRUNS {
            let input: Scalar = utils::random_scalar(true);
            let output = rc.compose(&rc.decompose(&input));

            assert_eq!(input, output);
        }
    }

    #[test]
    fn consistent_hash() {
        let rc = ReinforcedConcrete::new(&RC_BN_PARAMS);
        for _ in 0..TESTRUNS {
            let input1: Scalar = utils::random_scalar(true);
            let mut input2: Scalar;
            loop {
                input2 = utils::random_scalar(true);
                if input1 != input2 {
                    break;
                }
            }
            let input3: Scalar = utils::random_scalar(true);

            let h1 = rc.hash(&input1, &input3);
            let h2 = rc.hash(&input1, &input3);
            let h3 = rc.hash(&input2, &input3);
            assert_eq!(h1, h2);
            assert_ne!(h1, h3);
        }
    }

    #[test]
    fn kats() {
        let rc = ReinforcedConcrete::new(&RC_BN_PARAMS);
        let input: [Scalar; 3] = [Scalar::zero(), Scalar::one(), utils::from_u64(2)];
        let perm = rc.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x13f044e4afba717999129f20ff53f1fd0b3cf2601ce68048b503d8d37050d7b7").unwrap()
        );
        assert_eq!(
            perm[1],
            from_hex("0x20c6c94b8d312a2b78a5a3da7a922950c0a64c91340edd7db1a5fa954ad329fa").unwrap(),
        );
        assert_eq!(
            perm[2],
            from_hex("0x12f4d5c53ca385744684185fb2a26bbf54c93d1c2b7047cb891a5d6f7a9abc5b").unwrap(),
        );
    }
}

#[cfg(test)]
mod reinforced_concrete_tests_bls12 {
    use ff::{from_hex, Field};

    use crate::{
        fields::bls12::FpBLS12, reinforced_concrete::reinforced_concrete_instances::RC_BLS_PARAMS,
    };

    type Scalar = FpBLS12;

    use super::*;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let rc = ReinforcedConcrete::new(&RC_BLS_PARAMS);
        for _ in 0..TESTRUNS {
            let input1: [Scalar; 3] = [
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
            ];

            let mut input2: [Scalar; 3];
            loop {
                input2 = [
                    utils::random_scalar(true),
                    utils::random_scalar(true),
                    utils::random_scalar(true),
                ];
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = rc.permutation(&input1);
            let perm2 = rc.permutation(&input1);
            let perm3 = rc.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn compose() {
        let rc = ReinforcedConcrete::new(&RC_BLS_PARAMS);

        for _ in 0..TESTRUNS {
            let input: Scalar = utils::random_scalar(true);
            let output = rc.compose(&rc.decompose(&input));

            assert_eq!(input, output);
        }
    }

    #[test]
    fn consistent_hash() {
        let rc = ReinforcedConcrete::new(&RC_BLS_PARAMS);
        for _ in 0..TESTRUNS {
            let input1: Scalar = utils::random_scalar(true);
            let mut input2: Scalar;
            loop {
                input2 = utils::random_scalar(true);
                if input1 != input2 {
                    break;
                }
            }
            let input3: Scalar = utils::random_scalar(true);

            let h1 = rc.hash(&input1, &input3);
            let h2 = rc.hash(&input1, &input3);
            let h3 = rc.hash(&input2, &input3);
            assert_eq!(h1, h2);
            assert_ne!(h1, h3);
        }
    }

    #[test]
    fn kats() {
        let rc = ReinforcedConcrete::new(&RC_BLS_PARAMS);
        let input: [Scalar; 3] = [Scalar::zero(), Scalar::one(), utils::from_u64(2)];
        let perm = rc.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x0f628f13dc23302b0163eb586f9f939b483a7bd188e258b3b4631160eefd7a98").unwrap()
        );
        assert_eq!(
            perm[1],
            from_hex("0x50e932fb92013832944f5d36977f4900d4288bde1abe1f74569c8c71fbb04ba1").unwrap(),
        );
        assert_eq!(
            perm[2],
            from_hex("0x26c4a46a378a3226a677984ade7dc2f92139438f9fb47452e3bc64755591be62").unwrap(),
        );
    }
}
