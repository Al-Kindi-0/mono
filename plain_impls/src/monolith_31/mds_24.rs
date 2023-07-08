use crate::fields::f31::Field32;
use ff::PrimeField;

//row = [8, 59, 36, 24, 37, 15, 25, 27, 69, 95, 91, 20, 62, 64, 74, 52, 65, 53, 30, 81, 60, 4,
//            85, 29]
const MDS_FREQ_BLOCK_ONE: [i64; 6] = [125, 117, 278, 175, 239, 231];
const MDS_FREQ_BLOCK_TWO: [(i64, i64); 6] = [
    (-54, -5),
    (9, -38),
    (-6, -28),
    (-91, -28),
    (-9, -38),
    (54, -5),
];
const MDS_FREQ_BLOCK_THREE: [i64; 6] = [15, -19, 74, 23, 19, -15];
#[inline(always)]
#[allow(clippy::shadow_unrelated)]
pub(crate) fn mds_mult_24(state: [u64; 24]) -> [u64; 24] {
    let [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23] =
        state;

    let (u0, u1, u2) = fft4_real([s0, s6, s12, s18]);
    let (u4, u5, u6) = fft4_real([s1, s7, s13, s19]);
    let (u8, u9, u10) = fft4_real([s2, s8, s14, s20]);
    let (u12, u13, u14) = fft4_real([s3, s9, s15, s21]);
    let (u16, u17, u18) = fft4_real([s4, s10, s16, s22]);
    let (u20, u21, u22) = fft4_real([s5, s11, s17, s23]);

    let [v0, v4, v8, v12, v16, v20] = block1([u0, u4, u8, u12, u16, u20], MDS_FREQ_BLOCK_ONE);
    let [v1, v5, v9, v13, v17, v21] = block2([u1, u5, u9, u13, u17, u21], MDS_FREQ_BLOCK_TWO);
    let [v2, v6, v10, v14, v18, v22] = block3([u2, u6, u10, u14, u18, u22], MDS_FREQ_BLOCK_THREE);
    // The 4th block is not computed as it is similar to the 2nd one, up to complex conjugation,
    // and is, due to the use of the real FFT and iFFT, redundant.

    let [s0, s6, s12, s18] = ifft4_real((v0, v1, v2));
    let [s1, s7, s13, s19] = ifft4_real((v4, v5, v6));
    let [s2, s8, s14, s20] = ifft4_real((v8, v9, v10));
    let [s3, s9, s15, s21] = ifft4_real((v12, v13, v14));
    let [s4, s10, s16, s22] = ifft4_real((v16, v17, v18));
    let [s5, s11, s17, s23] = ifft4_real((v20, v21, v22));

    [
        s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19,
        s20, s21, s22, s23,
    ]
}

// We use the real FFT to avoid redundant computations. See https://www.mdpi.com/2076-3417/12/9/4700
#[inline(always)]
fn fft2_real(x: [u64; 2]) -> [i64; 2] {
    [(x[0] as i64 + x[1] as i64), (x[0] as i64 - x[1] as i64)]
}

#[inline(always)]
fn ifft2_real(y: [i64; 2]) -> [u64; 2] {
    [((y[0] + y[1]) / 2) as u64, ((y[0] - y[1]) / 2) as u64]
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
    let z0 = (y.0 + y.2) / 2;
    let z1 = (y.0 - y.2) / 2;

    let z2 = y.1 .0;
    let z3 = -y.1 .1;

    let [x0, x2] = ifft2_real([z0, z2]);
    let [x1, x3] = ifft2_real([z1, z3]);

    [x0, x1, x2, x3]
}

#[inline(always)]
fn block1(x: [i64; 6], y: [i64; 6]) -> [i64; 6] {
    let [x0, x1, x2, x3, x4, x5] = x;
    let [y0, y1, y2, y3, y4, y5] = y;

    let z0 = x0 * y0 + x1 * y5 + x2 * y4 + x3 * y3 + x4 * y2 + x5 * y1;
    let z1 = x0 * y1 + x1 * y0 + x2 * y5 + x3 * y4 + x4 * y3 + x5 * y2;
    let z2 = x0 * y2 + x1 * y1 + x2 * y0 + x3 * y5 + x4 * y4 + x5 * y3;
    let z3 = x0 * y3 + x1 * y2 + x2 * y1 + x3 * y0 + x4 * y5 + x5 * y4;
    let z4 = x0 * y4 + x1 * y3 + x2 * y2 + x3 * y1 + x4 * y0 + x5 * y5;
    let z5 = x0 * y5 + x1 * y4 + x2 * y3 + x3 * y2 + x4 * y1 + x5 * y0;

    [z0, z1, z2, z3, z4, z5]
}

#[inline(always)]
#[allow(clippy::shadow_unrelated)]
fn block2(x: [(i64, i64); 6], y: [(i64, i64); 6]) -> [(i64, i64); 6] {
    let [(x0r, x0i), (x1r, x1i), (x2r, x2i), (x3r, x3i), (x4r, x4i), (x5r, x5i)] = x;
    let [(y0r, y0i), (y1r, y1i), (y2r, y2i), (y3r, y3i), (y4r, y4i), (y5r, y5i)] = y;
    let x0s = x0r + x0i;
    let x1s = x1r + x1i;
    let x2s = x2r + x2i;
    let x3s = x3r + x3i;
    let x4s = x4r + x4i;
    let x5s = x5r + x5i;
    let y0s = y0r + y0i;
    let y1s = y1r + y1i;
    let y2s = y2r + y2i;
    let y3s = y3r + y3i;
    let y4s = y4r + y4i;
    let y5s = y5r + y5i;

    // Compute x0*y0 - I*x1*y5 - I*x2*y4 - I*x3*y3 - I*x4*y2 - I*x5*y1​ using Karatsuba
    let m0 = (x0r * y0r, x0i * y0i);
    let m1 = (x1r * y5r, x1i * y5i);
    let m2 = (x2r * y4r, x2i * y4i);
    let m3 = (x3r * y3r, x3i * y3i);
    let m4 = (x4r * y2r, x4i * y2i);
    let m5 = (x5r * y1r, x5i * y1i);

    let z0r = (m0.0 - m0.1)
        + (x1s * y5s - m1.0 - m1.1)
        + (x2s * y4s - m2.0 - m2.1)
        + (x3s * y3s - m3.0 - m3.1)
        + (x4s * y2s - m4.0 - m4.1)
        + (x5s * y1s - m5.0 - m5.1);
    let z0i = (x0s * y0s - m0.0 - m0.1)
        + (-m1.0 + m1.1)
        + (-m2.0 + m2.1)
        + (-m3.0 + m3.1)
        + (-m4.0 + m4.1)
        + (-m5.0 + m5.1);
    let z0 = (z0r, z0i);

    // x0*y1 + x1*y0 - I*x2*y5 - I*x3*y4 - I*x4*y3 - I*x5*y2
    let m0 = (x0r * y1r, x0i * y1i);
    let m1 = (x1r * y0r, x1i * y0i);
    let m2 = (x2r * y5r, x2i * y5i);
    let m3 = (x3r * y4r, x3i * y4i);
    let m4 = (x4r * y3r, x4i * y3i);
    let m5 = (x5r * y2r, x5i * y2i);
    let z1r = (m0.0 - m0.1)
        + (m1.0 - m1.1)
        + (x2s * y5s - m2.0 - m2.1)
        + (x3s * y4s - m3.0 - m3.1)
        + (x4s * y3s - m4.0 - m4.1)
        + (x5s * y2s - m5.0 - m5.1);
    let z1i = (x0s * y1s - m0.0 - m0.1)
        + (x1s * y0s - m1.0 - m1.1)
        + (-m2.0 + m2.1)
        + (-m3.0 + m3.1)
        + (-m4.0 + m4.1)
        + (-m5.0 + m5.1);
    let z1 = (z1r, z1i);

    // x0*y2 + x1*y1 + x2*y0 - I*x3*y5 - I*x4*y4 - I*x5*y3
    let m0 = (x0r * y2r, x0i * y2i);
    let m1 = (x1r * y1r, x1i * y1i);
    let m2 = (x2r * y0r, x2i * y0i);
    let m3 = (x3r * y5r, x3i * y5i);
    let m4 = (x4r * y4r, x4i * y4i);
    let m5 = (x5r * y3r, x5i * y3i);
    let z2r = (m0.0 - m0.1)
        + (m1.0 - m1.1)
        + (m2.0 - m2.1)
        + (x3s * y5s - m3.0 - m3.1)
        + (x4s * y4s - m4.0 - m4.1)
        + (x5s * y3s - m5.0 - m5.1);
    let z2i = (x0s * y2s - m0.0 - m0.1)
        + (x1s * y1s - m1.0 - m1.1)
        + (x2s * y0s - m2.0 - m2.1)
        + (-m3.0 + m3.1)
        + (-m4.0 + m4.1)
        + (-m5.0 + m5.1);
    let z2 = (z2r, z2i);

    // x0*y3 + x1*y2 + x2*y1 + x3*y0 - I*x4*y5 - I*x5*y4
    let m0 = (x0r * y3r, x0i * y3i);
    let m1 = (x1r * y2r, x1i * y2i);
    let m2 = (x2r * y1r, x2i * y1i);
    let m3 = (x3r * y0r, x3i * y0i);
    let m4 = (x4r * y5r, x4i * y5i);
    let m5 = (x5r * y4r, x5i * y4i);
    let z3r = (m0.0 - m0.1)
        + (m1.0 - m1.1)
        + (m2.0 - m2.1)
        + (m3.0 - m3.1)
        + (x4s * y5s - m4.0 - m4.1)
        + (x5s * y4s - m5.0 - m5.1);
    let z3i = (x0s * y3s - m0.0 - m0.1)
        + (x1s * y2s - m1.0 - m1.1)
        + (x2s * y1s - m2.0 - m2.1)
        + (x3s * y0s - m3.0 - m3.1)
        + (-m4.0 + m4.1)
        + (-m5.0 + m5.1);
    let z3 = (z3r, z3i);

    // x0*y4 + x1*y3 + x2*y2 + x3*y1 + x4*y0 - I*x5*y5
    let m0 = (x0r * y4r, x0i * y4i);
    let m1 = (x1r * y3r, x1i * y3i);
    let m2 = (x2r * y2r, x2i * y2i);
    let m3 = (x3r * y1r, x3i * y1i);
    let m4 = (x4r * y0r, x4i * y0i);
    let m5 = (x5r * y5r, x5i * y5i);
    let z4r = (m0.0 - m0.1)
        + (m1.0 - m1.1)
        + (m2.0 - m2.1)
        + (m3.0 - m3.1)
        + (m4.0 - m4.1)
        + (x5s * y5s - m5.0 - m5.1);
    let z4i = (x0s * y4s - m0.0 - m0.1)
        + (x1s * y3s - m1.0 - m1.1)
        + (x2s * y2s - m2.0 - m2.1)
        + (x3s * y1s - m3.0 - m3.1)
        + (x4s * y0s - m4.0 - m4.1)
        + (-m5.0 + m5.1);
    let z4 = (z4r, z4i);

    // x0*y5 + x1*y4 + x2*y3 + x3*y2 + x4*y1 + x5*y0
    let m0 = (x0r * y5r, x0i * y5i);
    let m1 = (x1r * y4r, x1i * y4i);
    let m2 = (x2r * y3r, x2i * y3i);
    let m3 = (x3r * y2r, x3i * y2i);
    let m4 = (x4r * y1r, x4i * y1i);
    let m5 = (x5r * y0r, x5i * y0i);
    let z5r = (m0.0 - m0.1)
        + (m1.0 - m1.1)
        + (m2.0 - m2.1)
        + (m3.0 - m3.1)
        + (m4.0 - m4.1)
        + (m5.0 - m5.1);
    let z5i = (x0s * y5s - m0.0 - m0.1)
        + (x1s * y4s - m1.0 - m1.1)
        + (x2s * y3s - m2.0 - m2.1)
        + (x3s * y2s - m3.0 - m3.1)
        + (x4s * y1s - m4.0 - m4.1)
        + (x5s * y0s - m5.0 - m5.1);
    let z5 = (z5r, z5i);

    [z0, z1, z2, z3, z4, z5]
}

#[inline(always)]
fn block3(x: [i64; 6], y: [i64; 6]) -> [i64; 6] {
    let [x0, x1, x2, x3, x4, x5] = x;
    let [y0, y1, y2, y3, y4, y5] = y;

    let z0 = x0 * y0 - x1 * y5 - x2 * y4 - x3 * y3 - x4 * y2 - x5 * y1;
    let z1 = x0 * y1 + x1 * y0 - x2 * y5 - x3 * y4 - x4 * y3 - x5 * y2;
    let z2 = x0 * y2 + x1 * y1 + x2 * y0 - x3 * y5 - x4 * y4 - x5 * y3;
    let z3 = x0 * y3 + x1 * y2 + x2 * y1 + x3 * y0 - x4 * y5 - x5 * y4;
    let z4 = x0 * y4 + x1 * y3 + x2 * y2 + x3 * y1 + x4 * y0 - x5 * y5;
    let z5 = x0 * y5 + x1 * y4 + x2 * y3 + x3 * y2 + x4 * y1 + x5 * y0;

    [z0, z1, z2, z3, z4, z5]
}

#[inline(always)]
pub fn mds_multiply<F: Field32 + PrimeField>(state: &mut [F; 24]) {
    mds_multiply_generic(state)
}

#[inline(always)]
pub fn mds_multiply_with_rc<F: Field32 + PrimeField>(
    state: &mut [F; 24],
    round_constants: &[F; 24],
) {
    mds_multiply_with_rc_generic(state, round_constants)
}

#[inline(always)]
pub fn mds_multiply_u64<F: Field32 + PrimeField>(state: &mut [u64; 24]) {
    mds_multiply_u64_generic::<F>(state)
}

#[inline(always)]
pub fn mds_multiply_with_rc_u64<F: Field32 + PrimeField>(
    state: &mut [u64; 24],
    round_constants: &[F; 24],
) {
    mds_multiply_with_rc_u64_generic(state, round_constants)
}

#[allow(unused)]
fn mds_multiply_generic<F: Field32 + PrimeField>(state: &mut [F; 24]) {
    let mut inner: [u64; 24] = [0; 24];
    for (i, b) in state.iter().enumerate() {
        inner[i] = b.to_u32() as u64;
    }

    inner = mds_mult_24(inner);

    for r in 0..24 {
        let s = inner[r];
        state[r] = F::from_u64(s);
    }
}

#[allow(unused)]
fn mds_multiply_with_rc_generic<F: Field32 + PrimeField>(
    state: &mut [F; 24],
    round_constants: &[F; 24],
) {
    let mut inner: [u64; 24] = [0; 24];
    for (i, b) in state.iter().enumerate() {
        inner[i] = b.to_u32() as u64;
    }

    inner = mds_mult_24(inner);

    for r in 0..16 {
        let s = inner[r] + round_constants[r].to_u32() as u64;
        state[r] = F::from_u64(s);
    }
}

#[allow(unused)]
fn mds_multiply_u64_generic<F: Field32 + PrimeField>(state: &mut [u64; 24]) {
    *state = mds_mult_24(*state);

    for el in state.iter_mut() {
        F::reduce64(el);
    }
}

#[allow(unused)]
fn mds_multiply_with_rc_u64_generic<F: Field32 + PrimeField>(
    state: &mut [u64; 24],
    round_constants: &[F; 24],
) {
    let lo = mds_mult_24(*state);

    for r in 0..24 {
        state[r] = lo[r] + round_constants[r].to_u32() as u64;
        F::reduce64(&mut state[r]);
    }
}
