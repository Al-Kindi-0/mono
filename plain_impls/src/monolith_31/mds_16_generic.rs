// from https://github.com/Neptune-Crypto/twenty-first/blob/master/twenty-first/src/shared_math/tip5.rs

use crate::fields::f31::Field32;
use ff::PrimeField;
use std::convert::TryInto;

/// The defining, first column of the (circulant) MDS matrix.
/// Derived from the SHA-256 hash of the ASCII string “Tip5” by dividing the digest into 16-bit
/// chunks.
pub const MDS_MATRIX_FIRST_COLUMN: [i64; 16] = [
    61402, 1108, 28750, 33823, 7454, 43244, 53865, 12034, 56951, 27521, 41351, 40901, 12021, 59689,
    26798, 17845,
];

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
    let mut lo: [u64; 16] = [0; 16];
    for (i, b) in state.iter().enumerate() {
        lo[i] = b.to_u32() as u64;
    }

    lo = fast_cyclomul16(&lo, MDS_MATRIX_FIRST_COLUMN);

    for r in 0..16 {
        let s = lo[r];
        state[r] = F::from_u64(s);
    }
}

#[allow(unused)]
fn mds_multiply_with_rc_generic<F: Field32 + PrimeField>(
    state: &mut [F; 16],
    round_constants: &[F; 16],
) {
    let mut lo: [u64; 16] = [0; 16];
    for (i, b) in state.iter().enumerate() {
        lo[i] = b.to_u32() as u64;
    }

    lo = fast_cyclomul16(&lo, MDS_MATRIX_FIRST_COLUMN);

    for r in 0..16 {
        let s = lo[r] + round_constants[r].to_u32() as u64;
        state[r] = F::from_u64(s);
    }
}

#[allow(unused)]
fn mds_multiply_u64_generic<F: Field32 + PrimeField>(state: &mut [u64; 16]) {
    *state = fast_cyclomul16(state, MDS_MATRIX_FIRST_COLUMN);

    for el in state.iter_mut() {
        F::reduce64(el);
    }
}

#[allow(unused)]
fn mds_multiply_with_rc_u64_generic<F: Field32 + PrimeField>(
    state: &mut [u64; 16],
    round_constants: &[F; 16],
) {
    let lo = fast_cyclomul16(state, MDS_MATRIX_FIRST_COLUMN);

    for r in 0..16 {
        state[r] = lo[r] + round_constants[r].to_u32() as u64;
        F::reduce64(&mut state[r]);
    }
}

#[inline(always)]
fn fast_cyclomul16(f: &[u64; 16], g: [i64; 16]) -> [u64; 16] {
    const N: usize = 8;
    let mut ff_lo = [0i64; N];
    let mut gg_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    let mut gg_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] as i64 + f[i + N] as i64;
        ff_hi[i] = f[i] as i64 - f[i + N] as i64;
        gg_lo[i] = g[i] + g[i + N];
        gg_hi[i] = g[i] - g[i + N];
    }

    let hh_lo = fast_cyclomul8(ff_lo, gg_lo);
    let hh_hi = complex_negacyclomul8(ff_hi, gg_hi);

    let mut hh = [0u64; 2 * N];
    for i in 0..N {
        hh[i] = ((hh_lo[i] + hh_hi[i]) >> 1) as u64;
        hh[i + N] = ((hh_lo[i] - hh_hi[i]) >> 1) as u64;
    }

    hh
}

#[inline(always)]
fn fast_cyclomul8(f: [i64; 8], g: [i64; 8]) -> [i64; 8] {
    const N: usize = 4;
    let mut ff_lo = [0i64; N];
    let mut gg_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    let mut gg_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
        gg_lo[i] = g[i] + g[i + N];
        gg_hi[i] = g[i] - g[i + N];
    }

    let hh_lo = fast_cyclomul4(ff_lo, gg_lo);
    let hh_hi = complex_negacyclomul4(ff_hi, gg_hi);

    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }

    hh
}

#[inline(always)]
fn fast_cyclomul4(f: [i64; 4], g: [i64; 4]) -> [i64; 4] {
    const N: usize = 2;
    let mut ff_lo = [0i64; N];
    let mut gg_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    let mut gg_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
        gg_lo[i] = g[i] + g[i + N];
        gg_hi[i] = g[i] - g[i + N];
    }

    let hh_lo = fast_cyclomul2(ff_lo, gg_lo);
    let hh_hi = complex_negacyclomul2(ff_hi, gg_hi);

    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }

    hh
}

#[inline(always)]
fn fast_cyclomul2(f: [i64; 2], g: [i64; 2]) -> [i64; 2] {
    let ff_lo = f[0] + f[1];
    let ff_hi = f[0] - f[1];
    let gg_lo = g[0] + g[1];
    let gg_hi = g[0] - g[1];

    let hh_lo = ff_lo * gg_lo;
    let hh_hi = ff_hi * gg_hi;

    let mut hh = [0i64; 2];
    hh[0] = (hh_lo + hh_hi) >> 1;
    hh[1] = (hh_lo - hh_hi) >> 1;

    hh
}

#[inline(always)]
fn complex_negacyclomul8(f: [i64; 8], g: [i64; 8]) -> [i64; 8] {
    const N: usize = 4;

    let mut f0 = [(0i64, 0i64); N];
    // let mut f1 = [(0i64,0i64); N];
    let mut g0 = [(0i64, 0i64); N];
    // let mut g1 = [(0i64,0i64); N];

    for i in 0..N {
        f0[i] = (f[i], -f[N + i]);
        // f1[i] = (f[i],  f[N+i]);
        g0[i] = (g[i], -g[N + i]);
        // g1[i] = (g[i],  g[N+i]);
    }

    let h0 = complex_karatsuba4(f0, g0);
    // h1 = complex_karatsuba(f1, g1)

    // h = a * h0 + b * h1
    // where a = 2^-1 * (i*X^(n/2) + 1)
    // and  b = 2^-1 * (-i*X^(n/2) + 1)

    let mut h = [0i64; 3 * N - 1];
    for i in 0..(2 * N - 1) {
        h[i] += h0[i].0;
        h[i + N] -= h0[i].1;
        // h[i] += h0[i].0 / 2
        // h[i+N] -= h0[i].1 / 2
        // h[i] += h1[i].0 / 2
        // h[i+N] -= h1[i].1 / 2
    }

    let mut hh = [0i64; 2 * N];
    for i in 0..(2 * N) {
        hh[i] += h[i];
    }
    for i in (2 * N)..(3 * N - 1) {
        hh[i - 2 * N] -= h[i];
    }

    hh
}

#[inline(always)]
fn complex_karatsuba4(f: [(i64, i64); 4], g: [(i64, i64); 4]) -> [(i64, i64); 7] {
    const N: usize = 2;

    let ff = complex_sum::<2>(f[..N].try_into().unwrap(), f[N..].try_into().unwrap());
    let gg = complex_sum::<2>(g[..N].try_into().unwrap(), g[N..].try_into().unwrap());

    let lo = complex_karatsuba2(f[..N].try_into().unwrap(), g[..N].try_into().unwrap());
    let hi = complex_karatsuba2(f[N..].try_into().unwrap(), g[N..].try_into().unwrap());

    let li = complex_diff::<3>(complex_karatsuba2(ff, gg), complex_sum::<3>(lo, hi));

    let mut result = [(0i64, 0i64); 4 * N - 1];
    for i in 0..(2 * N - 1) {
        result[i].0 = lo[i].0;
        result[i].1 = lo[i].1;
    }
    for i in 0..(2 * N - 1) {
        result[N + i].0 += li[i].0;
        result[N + i].1 += li[i].1;
    }
    for i in 0..(2 * N - 1) {
        result[2 * N + i].0 += hi[i].0;
        result[2 * N + i].1 += hi[i].1;
    }

    result
}

#[inline(always)]
fn complex_negacyclomul4(f: [i64; 4], g: [i64; 4]) -> [i64; 4] {
    const N: usize = 2;

    let mut f0 = [(0i64, 0i64); N];
    // let mut f1 = [(0i64,0i64); N];
    let mut g0 = [(0i64, 0i64); N];
    // let mut g1 = [(0i64,0i64); N];

    for i in 0..N {
        f0[i] = (f[i], -f[N + i]);
        // f1[i] = (f[i],  f[N+i]);
        g0[i] = (g[i], -g[N + i]);
        // g1[i] = (g[i],  g[N+i]);
    }

    let h0 = complex_karatsuba2(f0, g0);
    // h1 = complex_karatsuba(f1, g1)

    // h = a * h0 + b * h1
    // where a = 2^-1 * (i*X^(n/2) + 1)
    // and  b = 2^-1 * (-i*X^(n/2) + 1)

    let mut h = [0i64; 4 * N - 1];
    for i in 0..(2 * N - 1) {
        h[i] += h0[i].0;
        h[i + N] -= h0[i].1;
        // h[i] += h0[i].0 / 2
        // h[i+N] -= h0[i].1 / 2
        // h[i] += h1[i].0 / 2
        // h[i+N] -= h1[i].1 / 2
    }

    let mut hh = [0i64; 2 * N];
    for i in 0..(2 * N) {
        hh[i] += h[i];
    }
    for i in (2 * N)..(4 * N - 1) {
        hh[i - 2 * N] -= h[i];
    }

    hh
}

#[inline(always)]
fn complex_negacyclomul2(f: [i64; 2], g: [i64; 2]) -> [i64; 2] {
    let f0 = (f[0], -f[1]);
    let g0 = (g[0], -g[1]);

    let h0 = complex_product(f0, g0);

    [h0.0, -h0.1]
}

#[inline(always)]
fn complex_karatsuba2(f: [(i64, i64); 2], g: [(i64, i64); 2]) -> [(i64, i64); 3] {
    const N: usize = 1;

    let ff = (f[0].0 + f[1].0, f[0].1 + f[1].1);
    let gg = (g[0].0 + g[1].0, g[0].1 + g[1].1);

    let lo = complex_product(f[0], g[0]);
    let hi = complex_product(f[1], g[1]);

    let ff_times_gg = complex_product(ff, gg);
    let lo_plus_hi = (lo.0 + hi.0, lo.1 + hi.1);

    let li = (ff_times_gg.0 - lo_plus_hi.0, ff_times_gg.1 - lo_plus_hi.1);

    let mut result = [(0i64, 0i64); 4 * N - 1];
    result[0].0 += lo.0;
    result[0].1 += lo.1;
    result[N].0 += li.0;
    result[N].1 += li.1;
    result[2 * N].0 += hi.0;
    result[2 * N].1 += hi.1;

    result
}

#[inline(always)]
fn complex_sum<const N: usize>(f: [(i64, i64); N], g: [(i64, i64); N]) -> [(i64, i64); N] {
    let mut h = [(0i64, 0i64); N];
    for i in 0..N {
        h[i].0 = f[i].0 + g[i].0;
        h[i].1 = f[i].1 + g[i].1;
    }
    h
}

#[inline(always)]
fn complex_diff<const N: usize>(f: [(i64, i64); N], g: [(i64, i64); N]) -> [(i64, i64); N] {
    let mut h = [(0i64, 0i64); N];
    for i in 0..N {
        h[i].0 = f[i].0 - g[i].0;
        h[i].1 = f[i].1 - g[i].1;
    }
    h
}

#[inline(always)]
fn complex_product(f: (i64, i64), g: (i64, i64)) -> (i64, i64) {
    // don't karatsuba; this is faster
    (f.0 * g.0 - f.1 * g.1, f.0 * g.1 + f.1 * g.0)
}

///////////////////////////////////////////////////////////////////////////////
// test
///////////////////////////////////////////////////////////////////////////////

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
            61402, 17845, 26798, 59689, 12021, 40901, 41351, 27521, 56951, 12034, 53865, 43244,
            7454, 33823, 28750, 1108,
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

            // TODO adapt and check

            let output1 = matmul(&input, &mat);
            let mut output2 = input.to_owned();
            let mut output3 = input.to_owned();
            mds_multiply_with_rc_generic(&mut output2, &round_const);
            mds_multiply_generic(&mut output3);
            assert_eq!(output1, output2);
            assert_eq!(output1, output3);

            let mut output4 = [0u64; 16];
            for (src, des) in input.iter().zip(output4.iter_mut()) {
                *des = src.to_u32() as u64;
            }
            let mut output5 = output4.to_owned();
            mds_multiply_u64_generic::<Scalar>(&mut output4);
            mds_multiply_with_rc_u64_generic(&mut output5, &round_const);
            for (a, b) in output1.iter().zip(output4.iter()) {
                assert_eq!(a.to_u32(), *b as u32);
            }
            for (a, b) in output1.iter().zip(output5.iter()) {
                assert_eq!(a.to_u32(), *b as u32);
            }
        }
    }
}
