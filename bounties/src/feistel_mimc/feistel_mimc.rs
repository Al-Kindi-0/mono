use std::sync::Arc;

use ff::PrimeField;

use super::feistel_mimc_params::FeistelMimcParams;

#[derive(Clone, Debug)]
pub struct FeistelMimc<F: PrimeField> {
    pub(crate) params: Arc<FeistelMimcParams<F>>,
}

impl<F: PrimeField> FeistelMimc<F> {
    pub fn new(params: &Arc<FeistelMimcParams<F>>) -> Self {
        FeistelMimc {
            params: Arc::clone(params),
        }
    }

    fn sbox(&self, state_0: &F, round: usize) -> F {
        let mut input = *state_0;
        input.add_assign(&self.params.round_constants[round]);

        let mut input2 = input.clone();
        input2.square();
        match self.params.d {
            3 => {}
            5 => input2.square(),
            _ => assert!(false),
        }
        input2.mul_assign(&input);
        input2
    }

    fn round(&self, state: &mut [F; 2], round: usize) {
        let power = self.sbox(&state[0], round);
        state[1].add_assign(&power);
    }

    pub fn permutation(&self, input: &[F; 2]) -> [F; 2] {
        let mut current_state = input.to_owned();
        for r in 0..self.params.rounds - 1 {
            self.round(&mut current_state, r);
            current_state.swap(0, 1);
        }

        // finally without rotation
        self.round(&mut current_state, self.params.rounds - 1);

        current_state
    }
}

// #############################################################################

#[cfg(test)]
mod feistel_mimc_kats {
    use super::*;

    use crate::feistel_mimc::feistel_mimc_instances::*;
    use crate::fields::field::Fp;

    use ff::{from_hex, Field};

    type Scalar = Fp;

    #[test]
    fn easy1_kats() {
        let fm = FeistelMimc::new(&FM_PARAMS_EASY1);
        let input: [Scalar; 2] = [Scalar::zero(), Scalar::one()];
        let perm = fm.permutation(&input);
        assert_eq!(perm[0], from_hex("0x45b9c3d2a4bbdc9e").unwrap());
        assert_eq!(perm[1], from_hex("00ddd29700b3a6bc76").unwrap());
    }

    #[test]
    fn easy2_kats() {
        let fm = FeistelMimc::new(&FM_PARAMS_EASY2);
        let input: [Scalar; 2] = [Scalar::zero(), Scalar::one()];
        let perm = fm.permutation(&input);
        assert_eq!(perm[0], from_hex("0x788a3f6ea5a6a53e").unwrap());
        assert_eq!(perm[1], from_hex("0x25222a199c56f899").unwrap());
    }

    #[test]
    fn medium_kats() {
        let fm = FeistelMimc::new(&FM_PARAMS_MEDIUM);
        let input: [Scalar; 2] = [Scalar::zero(), Scalar::one()];
        let perm = fm.permutation(&input);
        assert_eq!(perm[0], from_hex("0x0bd74ec3122b04f6").unwrap());
        assert_eq!(perm[1], from_hex("0x67d7f4198480eaee").unwrap());
    }

    #[test]
    fn hard1_kats() {
        let fm = FeistelMimc::new(&FM_PARAMS_HARD1);
        let input: [Scalar; 2] = [Scalar::zero(), Scalar::one()];
        let perm = fm.permutation(&input);
        assert_eq!(perm[0], from_hex("0x299f5eaff8f41844").unwrap());
        assert_eq!(perm[1], from_hex("0xd7c9d66716f2cae7").unwrap());
    }

    #[test]
    fn hard2_kats() {
        let fm = FeistelMimc::new(&FM_PARAMS_HARD2);
        let input: [Scalar; 2] = [Scalar::zero(), Scalar::one()];
        let perm = fm.permutation(&input);
        assert_eq!(perm[0], from_hex("0xf874e35bbaf92376").unwrap());
        assert_eq!(perm[1], from_hex("0x72af6f65901ac3f1").unwrap());
    }
}
