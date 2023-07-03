use super::monolith_31::Monolith31;
use crate::fields::f31::Field32;
use ff::PrimeField;
use generic_array::typenum::{Unsigned, U8};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128, Shake128Reader,
};
use std::convert::TryInto;

// weird workaround for not allowed const generics in arrays
pub(crate) type BarsType = U8; // U8 = 8 bars per layer

#[derive(Clone, Debug)]
pub struct Monolith31Params<F: Field32 + PrimeField, const T: usize> {
    pub(crate) round_constants: Vec<[F; T]>,
    #[allow(unused)]
    pub(crate) mds: [[F; T]; T],
    pub(crate) lookup1: Vec<u16>,
    pub(crate) lookup2: Vec<u16>,
}

impl<F: Field32 + PrimeField, const T: usize> Monolith31Params<F, T> {
    pub const R: usize = 6;
    pub const BARS: usize = BarsType::USIZE; // Change BarsType to change the value
    pub const INIT_SHAKE: &'static str = "Monolith";

    pub fn new() -> Self {
        assert_eq!(F::NUM_BITS, 31);
        assert!(T >= 16);
        assert_eq!(T % 4, 0);
        let round_constants = Self::instantiate_rc();
        let lookup1 = Self::instantiate_lookup1();
        let lookup2 = Self::instantiate_lookup2();
        let mds = Self::get_mds();

        Monolith31Params {
            round_constants,
            mds,
            lookup1,
            lookup2,
        }
    }

    fn bar0_8(limb: u8) -> u8 {
        let limbl1 = (limb >> 7) | (limb << 1); //left rot by 1
        let limbl2 = (limb >> 6) | (limb << 2); //left rot by 2
        let limbl3 = (limb >> 5) | (limb << 3); //left rot by 3

        //yi = xi +  (1 + x{i+1}) * x{i+2} * x{i+3}
        let tmp = limb ^ !limbl1 & limbl2 & limbl3;
        (tmp >> 7) | (tmp << 1) // Final rotate for less S(x) = x
    }

    // We keep a u16 lookup here
    fn instantiate_lookup1() -> Vec<u16> {
        (0..=u16::MAX)
            .map(|i| {
                let a = (i >> 8) as u8;
                let b = i as u8;
                ((Self::bar0_8(a) as u16) << 8) | Self::bar0_8(b) as u16
            })
            .collect()
    }

    // We keep a u16 lookup here
    fn instantiate_lookup2() -> Vec<u16> {
        (0..(1 << 15))
            .map(|i| {
                let a = (i >> 8) as u8;
                let b = i as u8;
                ((Monolith31::<F, T>::bar1_7(a) as u16) << 8) | Self::bar0_8(b) as u16
            })
            .collect()
    }

    fn random_field_element(shake: &mut Shake128Reader) -> F {
        loop {
            let mut rnd = [0u8; 4];
            shake.read(&mut rnd);
            let repr = F::Repr::from(u32::from_le_bytes(rnd) as u64);
            let res = F::from_repr(repr);
            if let Ok(r) = res {
                return r;
            }
        }
    }

    fn circ_mat(row: &[u64; T]) -> [[F; T]; T] {
        let mut mat = [[F::zero(); T]; T];
        let mut rot: Vec<F> = row.iter().map(|i| F::from_u64(*i)).collect();
        mat[0].copy_from_slice(&rot);
        for row in mat.iter_mut().skip(1) {
            rot.rotate_right(1);
            row.copy_from_slice(&rot);
        }
        mat
    }

    fn get_mds() -> [[F; T]; T] {
        if T == 16 {
            let row = [
                61402, 17845, 26798, 59689, 12021, 40901, 41351, 27521, 56951, 12034, 53865, 43244,
                7454, 33823, 28750, 1108,
            ];
            Self::circ_mat(row.as_ref().try_into().unwrap())
        } else {
            let mut shake = Shake128::default();
            shake.update(Self::INIT_SHAKE.as_bytes());
            shake.update(&[T as u8, Self::R as u8]);
            shake.update(&u32::to_le_bytes(F::char().as_ref()[0] as u32));
            shake.update(&[16, 15]);
            shake.update("MDS".as_bytes());
            let mut shake = shake.finalize_xof();
            Self::cauchy_mds_matrix(&mut shake)
        }
    }

    //----------------------------------------------------------------
    // construct cauchy matrix:
    // get t distinct x_i
    // get t distinct y_i
    // if for all i, j: x_i + y_i != 0
    // then: a_i_j = (x_i + y_j)^-1 is MDS
    // construct:
    // - sample s+t y_i with floor(log_2(p)) - 1 bits
    // - set x_i to be the least significant ceil(log_2(p)) - r_bits - 2 bits of x_i
    // - x_i need to be distinct -> y_i also distinct
    // - requires sampling of s+t random values of size floor(log_2(p))
    fn cauchy_mds_matrix(shake: &mut Shake128Reader) -> [[F; T]; T] {
        let mut p = F::char().as_ref()[0].to_owned() as u32;
        let mut tmp = 0u32;
        while p != 0 {
            tmp += 1;
            p >>= 1;
        }
        // reduce two: to prevent that x_i + y_i = 0 for MDS
        let x_mask = (1 << (tmp - 7 - 2)) - 1;
        let y_mask = ((1 << tmp) - 1) >> 2;

        let mut res = [[F::zero(); T]; T];

        let y = Self::get_random_yi(shake, x_mask, y_mask);
        let mut x = y.to_owned();
        x.iter_mut().for_each(|xi| *xi &= x_mask);

        for (i, xi) in x.iter().enumerate() {
            for (j, yj) in y.iter().enumerate() {
                res[i][j] = F::from_u32(xi + yj).inverse().unwrap();
            }
        }

        res
    }

    fn get_random_yi(shake: &mut Shake128Reader, x_mask: u32, y_mask: u32) -> [u32; T] {
        let mut res = [0; T];
        for i in 0..T {
            let mut valid = false;
            while !valid {
                let mut rand = [0u8; 4];
                shake.read(&mut rand);
                let y_i = u32::from_be_bytes(rand) & y_mask;
                // check distinct x_i
                let x_i = y_i & x_mask;
                valid = true;
                for r in res.iter().take(i) {
                    if r & x_mask == x_i {
                        valid = false;
                        break;
                    }
                }
                if valid {
                    res[i] = y_i;
                }
            }
        }

        res
    }

    fn init_shake() -> Shake128Reader {
        // RCp<t><r><p in little endian><partition of sboxes>
        let mut shake = Shake128::default();
        shake.update(Self::INIT_SHAKE.as_bytes());
        shake.update(&[T as u8, Self::R as u8]);
        shake.update(&u32::to_le_bytes(F::char().as_ref()[0] as u32));
        shake.update(&[8, 8, 8, 7]);
        shake.finalize_xof()
    }

    fn instantiate_rc() -> Vec<[F; T]> {
        let mut shake = Self::init_shake();
        (0..Self::R - 1)
            .map(|_| {
                let mut rc = [F::zero(); T];
                rc.iter_mut()
                    .for_each(|el| *el = Self::random_field_element(&mut shake));
                rc
            })
            .collect()
    }

    pub fn get_t(&self) -> usize {
        T
    }

    pub fn get_rounds(&self) -> usize {
        Self::R
    }

    pub fn get_rc(&self, i: usize) -> [F; T] {
        self.round_constants[i]
    }
}

impl<F: PrimeField + Field32, const T: usize> Default for Monolith31Params<F, T> {
    fn default() -> Self {
        Self::new()
    }
}
