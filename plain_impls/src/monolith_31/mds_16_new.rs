use crate::fields::f31::Field32;
use ff::PrimeField;

const MDS_FREQ_BLOCK_ONE: [i64; 4] = [32, 24, 8, 32];
const MDS_FREQ_BLOCK_TWO: [(i64, i64); 4] = [(8, -16), (32, -14), (-1, 8), (16, 1)];
const MDS_FREQ_BLOCK_THREE: [i64; 4] = [-2, 9, -2, 4];

#[inline(always)]
#[allow(clippy::shadow_unrelated)]
pub(crate) fn mds_multiply_freq(state: [u64; 16]) -> [u64; 16] {
    let [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15] = state;

    let (u0, u1, u2) = fft4_real([s0, s4, s8, s12]);
    let (u4, u5, u6) = fft4_real([s1, s5, s9, s13]);
    let (u8, u9, u10) = fft4_real([s2, s6, s10, s14]);
    let (u12, u13, u14) = fft4_real([s3, s7, s11, s15]);

    let [v0, v4, v8, v12] = block1([u0, u4, u8, u12], MDS_FREQ_BLOCK_ONE);
    let [v1, v5, v9, v13] = block2([u1, u5, u9, u13], MDS_FREQ_BLOCK_TWO);
    let [v2, v6, v10, v14] = block3([u2, u6, u10, u14], MDS_FREQ_BLOCK_THREE);
    // The 4th block is not computed as it is similar to the 2nd one, up to complex conjugation,
    // and is, due to the use of the real FFT and iFFT, redundant.

    let [s0, s4, s8, s12] = ifft4_real((v0, v1, v2));
    let [s1, s5, s9, s13] = ifft4_real((v4, v5, v6));
    let [s2, s6, s10, s14] = ifft4_real((v8, v9, v10));
    let [s3, s7, s11, s15] = ifft4_real((v12, v13, v14));

    [
        s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15,
    ]
}

// We use the real FFT to avoid redundant computations. See https://www.mdpi.com/2076-3417/12/9/4700
#[inline(always)]
fn fft2_real(x: [u64; 2]) -> [i64; 2] {
    [(x[0] as i64 + x[1] as i64), (x[0] as i64 - x[1] as i64)]
}

#[inline(always)]
fn ifft2_real(y: [i64; 2]) -> [u64; 2] {
    // We avoid divisions by 2 by appropriately scaling the MDS matrix constants.
    [(y[0] + y[1]) as u64, (y[0] - y[1]) as u64]
}

#[inline(always)]
fn fft4_real(x: [u64; 4]) -> (i64, (i64, i64), i64) {
    let [z0, z2] = fft2_real([x[0], x[2]]);
    let [z1, z3] = fft2_real([x[1], x[3]]);
    let y0 = z0 + z1;
    let y1 = (z2, -z3);
    let y2 = z0 - z1;
    (y0, y1, y2)
}

#[inline(always)]
fn ifft4_real(y: (i64, (i64, i64), i64)) -> [u64; 4] {
    // In calculating 'z0' and 'z1', division by 2 is avoided by appropriately scaling
    // the MDS matrix constants.
    let z0 = y.0 + y.2;
    let z1 = y.0 - y.2;
    let z2 = y.1 .0;
    let z3 = -y.1 .1;

    let [x0, x2] = ifft2_real([z0, z2]);
    let [x1, x3] = ifft2_real([z1, z3]);

    [x0, x1, x2, x3]
}

#[inline(always)]
fn block1(x: [i64; 4], y: [i64; 4]) -> [i64; 4] {
    let [x0, x1, x2, x3] = x;
    let [y0, y1, y2, y3] = y;
    let z0 = x0 * y0 + x1 * y3 + x2 * y2 + x3 * y1;
    let z1 = x0 * y1 + x1 * y0 + x2 * y3 + x3 * y2;
    let z2 = x0 * y2 + x1 * y1 + x2 * y0 + x3 * y3;
    let z3 = x0 * y3 + x1 * y2 + x2 * y1 + x3 * y0;

    [z0, z1, z2, z3]
}

#[inline(always)]
#[allow(clippy::shadow_unrelated)]
fn block2(x: [(i64, i64); 4], y: [(i64, i64); 4]) -> [(i64, i64); 4] {
    let [(x0r, x0i), (x1r, x1i), (x2r, x2i), (x3r, x3i)] = x;
    let [(y0r, y0i), (y1r, y1i), (y2r, y2i), (y3r, y3i)] = y;
    let x0s = x0r + x0i;
    let x1s = x1r + x1i;
    let x2s = x2r + x2i;
    let x3s = x3r + x3i;
    let y0s = y0r + y0i;
    let y1s = y1r + y1i;
    let y2s = y2r + y2i;
    let y3s = y3r + y3i;

    // Compute x0​y0​−ix1​y3​−ix2​y2​−ix3​y1​ using Karatsuba
    let m0 = (x0r * y0r, x0i * y0i);
    let m1 = (x1r * y3r, x1i * y3i);
    let m2 = (x2r * y2r, x2i * y2i);
    let m3 = (x3r * y1r, x3i * y1i);
    let z0r = (m0.0 - m0.1)
        + (x1s * y3s - m1.0 - m1.1)
        + (x2s * y2s - m2.0 - m2.1)
        + (x3s * y1s - m3.0 - m3.1);
    let z0i = (x0s * y0s - m0.0 - m0.1) + (-m1.0 + m1.1) + (-m2.0 + m2.1) + (-m3.0 + m3.1);
    let z0 = (z0r, z0i);

    // Compute x0​y1​+x1​y0​−ix2​y3​−ix3​y2​ using Karatsuba
    let m0 = (x0r * y1r, x0i * y1i);
    let m1 = (x1r * y0r, x1i * y0i);
    let m2 = (x2r * y3r, x2i * y3i);
    let m3 = (x3r * y2r, x3i * y2i);
    let z1r = (m0.0 - m0.1) + (m1.0 - m1.1) + (x2s * y3s - m2.0 - m2.1) + (x3s * y2s - m3.0 - m3.1);
    let z1i =
        (x0s * y1s - m0.0 - m0.1) + (x1s * y0s - m1.0 - m1.1) + (-m2.0 + m2.1) + (-m3.0 + m3.1);
    let z1 = (z1r, z1i);

    // Compute x0​y2​+x1​y1​+x2​y0​−ix3​y3​​ using Karatsuba
    let m0 = (x0r * y2r, x0i * y2i);
    let m1 = (x1r * y1r, x1i * y1i);
    let m2 = (x2r * y0r, x2i * y0i);
    let m3 = (x3r * y3r, x3i * y3i);
    let z2r = (m0.0 - m0.1) + (m1.0 - m1.1) + (m2.0 - m2.1) + (x3s * y3s - m3.0 - m3.1);
    let z2i = (x0s * y2s - m0.0 - m0.1)
        + (x1s * y1s - m1.0 - m1.1)
        + (x2s * y0s - m2.0 - m2.1)
        + (-m3.0 + m3.1);
    let z2 = (z2r, z2i);

    // Compute x0​y3​+x1​y2​+x2​y1​+x3​y0​​​ using Karatsuba
    let m0 = (x0r * y3r, x0i * y3i);
    let m1 = (x1r * y2r, x1i * y2i);
    let m2 = (x2r * y1r, x2i * y1i);
    let m3 = (x3r * y0r, x3i * y0i);
    let z3r = (m0.0 - m0.1) + (m1.0 - m1.1) + (m2.0 - m2.1) + (m3.0 - m3.1);
    let z3i = (x0s * y3s - m0.0 - m0.1)
        + (x1s * y2s - m1.0 - m1.1)
        + (x2s * y1s - m2.0 - m2.1)
        + (x3s * y0s - m3.0 - m3.1);
    let z3 = (z3r, z3i);

    [z0, z1, z2, z3]
}

#[inline(always)]
fn block3(x: [i64; 4], y: [i64; 4]) -> [i64; 4] {
    let [x0, x1, x2, x3] = x;
    let [y0, y1, y2, y3] = y;

    let z0 = x0 * y0 - x1 * y3 - x2 * y2 - x3 * y1;
    let z1 = x0 * y1 + x1 * y0 - x2 * y3 - x3 * y2;
    let z2 = x0 * y2 + x1 * y1 + x2 * y0 - x3 * y3;
    let z3 = x0 * y3 + x1 * y2 + x2 * y1 + x3 * y0;

    [z0, z1, z2, z3]
}

#[inline(always)]
pub fn mds_multiply<F: Field32 + PrimeField>(state: &mut [F; 16]) {
    mds_multiply_generic(state)
}

#[inline(always)]
pub fn mds_multiply_with_rc<F: Field32 + PrimeField>(
    state: &mut [F; 16],
    round_constants: &[F; 16],
) {
    mds_multiply_with_rc_generic(state, round_constants)
}

#[inline(always)]
pub fn mds_multiply_u64<F: Field32 + PrimeField>(state: &mut [u64; 16]) {
    mds_multiply_u64_generic::<F>(state)
}

#[inline(always)]
pub fn mds_multiply_with_rc_u64<F: Field32 + PrimeField>(
    state: &mut [u64; 16],
    round_constants: &[F; 16],
) {
    mds_multiply_with_rc_u64_generic(state, round_constants)
}

#[allow(unused)]
fn mds_multiply_generic<F: Field32 + PrimeField>(state: &mut [F; 16]) {
    let mut inner: [u64; 16] = [0; 16];
    for (i, b) in state.iter().enumerate() {
        inner[i] = b.to_u32() as u64;
    }

    inner = mds_multiply_freq(inner);

    for r in 0..16 {
        let s = inner[r];
        state[r] = F::from_u64(s);
    }
}

#[allow(unused)]
fn mds_multiply_with_rc_generic<F: Field32 + PrimeField>(
    state: &mut [F; 16],
    round_constants: &[F; 16],
) {
    let mut inner: [u64; 16] = [0; 16];
    for (i, b) in state.iter().enumerate() {
        inner[i] = b.to_u32() as u64;
    }

    inner = mds_multiply_freq(inner);

    for r in 0..16 {
        let s = inner[r] + round_constants[r].to_u32() as u64;
        state[r] = F::from_u64(s);
    }
}

#[allow(unused)]
fn mds_multiply_u64_generic<F: Field32 + PrimeField>(state: &mut [u64; 16]) {
    *state = mds_multiply_freq(*state);

    for el in state.iter_mut() {
        F::reduce64(el);
    }
}

#[allow(unused)]
fn mds_multiply_with_rc_u64_generic<F: Field32 + PrimeField>(
    state: &mut [u64; 16],
    round_constants: &[F; 16],
) {
    let lo = mds_multiply_freq(*state);

    for r in 0..16 {
        state[r] = lo[r] + round_constants[r].to_u32() as u64;
        F::reduce64(&mut state[r]);
    }
}

#[cfg(test)]
mod mds_tests {
    use super::*;
    use crate::fields::{f31::F31, utils};
    use ff::Field;

    static TESTRUNS: usize = 5;
    type Scalar = F31;

    fn matmul(input: &[Scalar], mat: &[Vec<Scalar>]) -> Vec<Scalar> {
        let t = mat.len();
        debug_assert!(t == input.len());
        let mut out = vec![Scalar::zero(); t];
        for row in 0..t {
            for (col, inp) in input.iter().enumerate().take(t) {
                let mut tmp = mat[row][col];
                tmp.mul_assign(inp);
                out[row].add_assign(&tmp);
            }
        }
        out
    }

    fn circ_mat(row: &[u32]) -> Vec<Vec<Scalar>> {
        let t = row.len();
        let mut mat: Vec<Vec<Scalar>> = Vec::with_capacity(t);
        let mut rot: Vec<Scalar> = row.iter().map(|i| Scalar::from_u32(*i)).collect();
        mat.push(rot.clone());
        for _ in 1..t {
            rot.rotate_right(1);
            mat.push(rot.clone());
        }
        mat
    }

    #[test]
    fn kats() {
        let row = [
            38, 29, 18, 1, 18, 20, 7, 1, 22, 27, 2, 29, 50, 52, 5, 65
        ];
        let mat = circ_mat(&row);
        let round_const = [Scalar::zero(); 16];

        for _ in 0..TESTRUNS {
            let input: [Scalar; 16] = [
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
                utils::random_scalar(true),
            ];


            let output1 = matmul(&input, &mat);
            let mut output2 = input.to_owned();
            let mut output4 = input.to_owned();
            mds_multiply_with_rc_generic(&mut output2, &round_const);
            mds_multiply_generic(&mut output4);
            assert_eq!(output1, output2);
            assert_eq!(output1, output4);

            let mut output6 = [0u64; 16];
            for (src, des) in input.iter().zip(output6.iter_mut()) {
                *des = src.to_u32() as u64;
            }
            let mut output8 = output6.to_owned();
            mds_multiply_with_rc_u64_generic(&mut output6, &round_const);
            mds_multiply_u64_generic::<Scalar>(&mut output8);
            for (a, b) in output1.iter().zip(output6.iter()) {
                assert_eq!(a.to_u32(), *b as u32);
            }
            for (a, b) in output1.iter().zip(output8.iter()) {
                assert_eq!(a.to_u32(), *b as u32);
            }
        }
    }
}
