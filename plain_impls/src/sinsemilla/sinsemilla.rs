//! The Sinsemilla hash function and commitment scheme.

use super::sinsemilla_s::SINSEMILLA_S;
use halo2::arithmetic::CurveExt;
use pasta_curves::{pallas, Ep};
use subtle::CtOption;

use super::spec::extract_p_bottom;
use super::util::gen_const_array;

use super::addition::IncompletePoint;

use super::constants::*;

/// $\ell^\mathsf{Orchard}_\mathsf{Merkle}$
pub const L_ORCHARD_MERKLE: usize = 255;

/// SWU hash-to-curve personalization for the note commitment generator
pub const NOTE_COMMITMENT_PERSONALIZATION: &str = "z.cash:Orchard-NoteCommit";

/// SWU hash-to-curve personalization for the IVK commitment generator
pub const COMMIT_IVK_PERSONALIZATION: &str = "z.cash:Orchard-CommitIvk";

/// SWU hash-to-curve personalization for the Merkle CRH generator
pub const MERKLE_CRH_PERSONALIZATION: &str = "z.cash:Orchard-MerkleCRH";

pub(crate) fn lebs2ip_k(bits: &[bool]) -> u32 {
    assert!(bits.len() == K);
    bits.iter()
        .enumerate()
        .fold(0u32, |acc, (i, b)| acc + if *b { 1 << i } else { 0 })
}

/// The sequence of K bits in little-endian order representing an integer
/// up to `2^K` - 1.
pub fn i2lebsp_k(int: usize) -> [bool; K] {
    assert!(int < (1 << K));
    gen_const_array(|mask: usize| (int & (1 << mask)) != 0)
}

/// Pads the given iterator (which MUST have length $\leq K * C$) with zero-bits to a
/// multiple of $K$ bits.
struct Pad<I: Iterator<Item = bool>> {
    /// The iterator we are padding.
    inner: I,
    /// The measured length of the inner iterator.
    ///
    /// This starts as a lower bound, and will be accurate once `padding_left.is_some()`.
    len: usize,
    /// The amount of padding that remains to be emitted.
    padding_left: Option<usize>,
}

impl<I: Iterator<Item = bool>> Pad<I> {
    fn new(inner: I) -> Self {
        Pad {
            inner,
            len: 0,
            padding_left: None,
        }
    }
}

impl<I: Iterator<Item = bool>> Iterator for Pad<I> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we have identified the required padding, the inner iterator has ended,
            // and we will never poll it again.
            if let Some(n) = self.padding_left.as_mut() {
                if *n == 0 {
                    // Either we already emitted all necessary padding, or there was no
                    // padding required.
                    break None;
                } else {
                    // Emit the next padding bit.
                    *n -= 1;
                    break Some(false);
                }
            } else if let Some(ret) = self.inner.next() {
                // We haven't reached the end of the inner iterator yet.
                self.len += 1;
                assert!(self.len <= K * C);
                break Some(ret);
            } else {
                // Inner iterator just ended, so we now know its length.
                let rem = self.len % K;
                if rem > 0 {
                    // The inner iterator requires padding in the range [1,K).
                    self.padding_left = Some(K - rem);
                } else {
                    // No padding required.
                    self.padding_left = Some(0);
                }
            }
        }
    }
}

/// A domain in which $\mathsf{SinsemillaHashToPoint}$ and $\mathsf{SinsemillaHash}$ can
/// be used.
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct HashDomain {
    pub(crate) Q: pallas::Point,
}

impl HashDomain {
    /// Constructs a new `HashDomain` with a specific prefix string.
    pub fn new(domain: &str) -> Self {
        HashDomain {
            Q: pallas::Point::hash_to_curve(Q_PERSONALIZATION)(domain.as_bytes()),
        }
    }

    /// $\mathsf{SinsemillaHashToPoint}$ from [§ 5.4.1.9][concretesinsemillahash].
    ///
    /// [concretesinsemillahash]: https://zips.z.cash/protocol/nu5.pdf#concretesinsemillahash
    pub(crate) fn hash_to_point(&self, msg: impl Iterator<Item = bool>) -> CtOption<pallas::Point> {
        // self.hash_to_point_inner(msg).into()
        self.hash_to_point_inner_lookup(msg).into() // I added the lookup
    }

    #[allow(non_snake_case)]
    #[cfg(test)]
    fn hash_to_point_inner(&self, msg: impl Iterator<Item = bool>) -> IncompletePoint {
        let padded: Vec<_> = Pad::new(msg).collect();

        let hasher_S = pallas::Point::hash_to_curve(S_PERSONALIZATION);
        let S = |chunk: &[bool]| hasher_S(&lebs2ip_k(chunk).to_le_bytes());

        padded
            .chunks(K)
            .fold(IncompletePoint::from(self.Q), |acc, chunk| {
                (acc + S(chunk)) + acc
            })
    }

    // implemented by me
    #[allow(non_snake_case)]
    fn hash_to_point_inner_lookup(&self, msg: impl Iterator<Item = bool>) -> IncompletePoint {
        let padded: Vec<_> = Pad::new(msg).collect();

        padded
            .chunks(K)
            .fold(IncompletePoint::from(self.Q), |acc, chunk| {
                let lookedup = SINSEMILLA_S[lebs2ip_k(chunk) as usize];
                let point = Ep::new_jacobian(lookedup.0, lookedup.1, pallas::Base::one()).unwrap();
                (acc + point) + acc
            })
    }

    /// $\mathsf{SinsemillaHash}$ from [§ 5.4.1.9][concretesinsemillahash].
    ///
    /// [concretesinsemillahash]: https://zips.z.cash/protocol/nu5.pdf#concretesinsemillahash
    pub fn hash(&self, msg: impl Iterator<Item = bool>) -> CtOption<pallas::Base> {
        extract_p_bottom(self.hash_to_point(msg))
    }

    // /// Returns the Sinsemilla $Q$ constant for this domain.
    // #[cfg(test)]
    // #[allow(non_snake_case)]
    // pub(crate) fn Q(&self) -> pallas::Point {
    //     self.Q
    // }
}

/// A domain in which $\mathsf{SinsemillaCommit}$ and $\mathsf{SinsemillaShortCommit}$ can
/// be used.
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct CommitDomain {
    pub(crate) M: HashDomain,
    R: pallas::Point,
}

impl CommitDomain {
    /// Constructs a new `CommitDomain` with a specific prefix string.
    pub fn new(domain: &str) -> Self {
        let m_prefix = format!("{}-M", domain);
        let r_prefix = format!("{}-r", domain);
        let hasher_r = pallas::Point::hash_to_curve(&r_prefix);
        CommitDomain {
            M: HashDomain::new(&m_prefix),
            R: hasher_r(&[]),
        }
    }

    /// $\mathsf{SinsemillaCommit}$ from [§ 5.4.8.4][concretesinsemillacommit].
    ///
    /// [concretesinsemillacommit]: https://zips.z.cash/protocol/nu5.pdf#concretesinsemillacommit
    #[allow(non_snake_case)]
    pub(crate) fn commit(
        &self,
        msg: impl Iterator<Item = bool>,
        r: &pallas::Scalar,
    ) -> CtOption<pallas::Point> {
        // (self.M.hash_to_point_inner(msg) + self.R * r).into()
        (self.M.hash_to_point_inner_lookup(msg) + self.R * r).into()
    }

    /// $\mathsf{SinsemillaShortCommit}$ from [§ 5.4.8.4][concretesinsemillacommit].
    ///
    /// [concretesinsemillacommit]: https://zips.z.cash/protocol/nu5.pdf#concretesinsemillacommit
    pub fn short_commit(
        &self,
        msg: impl Iterator<Item = bool>,
        r: &pallas::Scalar,
    ) -> CtOption<pallas::Base> {
        extract_p_bottom(self.commit(msg, r))
    }

    // /// Returns the Sinsemilla $R$ constant for this domain.
    // #[cfg(test)]
    // #[allow(non_snake_case)]
    // pub(crate) fn R(&self) -> pallas::Point {
    //     self.R
    // }
}

#[cfg(test)]
mod tests {
    use super::{i2lebsp_k, lebs2ip_k, Pad, K};
    use random::{self, rngs::OsRng, Rng};

    #[test]
    fn pad() {
        assert_eq!(Pad::new([].iter().cloned()).collect::<Vec<_>>(), vec![]);
        assert_eq!(
            Pad::new([true].iter().cloned()).collect::<Vec<_>>(),
            vec![true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Pad::new([true, true].iter().cloned()).collect::<Vec<_>>(),
            vec![true, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Pad::new([true, true, true].iter().cloned()).collect::<Vec<_>>(),
            vec![true, true, true, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Pad::new(
                [true, true, false, true, false, true, false, true, false, true]
                    .iter()
                    .cloned()
            )
            .collect::<Vec<_>>(),
            vec![true, true, false, true, false, true, false, true, false, true]
        );
        assert_eq!(
            Pad::new(
                [true, true, false, true, false, true, false, true, false, true, true]
                    .iter()
                    .cloned()
            )
            .collect::<Vec<_>>(),
            vec![
                true, true, false, true, false, true, false, true, false, true, true, false, false,
                false, false, false, false, false, false, false
            ]
        );
    }

    #[test]
    fn lebs2ip_k_round_trip() {
        let mut rng = OsRng;
        {
            let int = rng.gen_range(0..(1 << K));
            assert_eq!(lebs2ip_k(&i2lebsp_k(int)) as usize, int);
        }

        assert_eq!(lebs2ip_k(&i2lebsp_k(0)) as usize, 0);
        assert_eq!(lebs2ip_k(&i2lebsp_k((1 << K) - 1)) as usize, (1 << K) - 1);
    }

    #[test]
    fn i2lebsp_k_round_trip() {
        {
            let bitstring = (0..K).map(|_| rand::random()).collect::<Vec<_>>();
            assert_eq!(
                i2lebsp_k(lebs2ip_k(&bitstring) as usize).to_vec(),
                bitstring
            );
        }

        {
            let bitstring = [false; K];
            assert_eq!(
                i2lebsp_k(lebs2ip_k(&bitstring) as usize).to_vec(),
                bitstring
            );
        }

        {
            let bitstring = [true; K];
            assert_eq!(
                i2lebsp_k(lebs2ip_k(&bitstring) as usize).to_vec(),
                bitstring
            );
        }
    }
}

// implemented by me
#[cfg(test)]
mod sinsemilla_tests {
    use super::*;
    use group::ff::Field;
    use group::ff::PrimeFieldBits;
    use pasta_curves::pallas::Base;
    use random::thread_rng;
    use std::iter;

    static TESTRUNS: usize = 5;

    #[test]
    fn hash_to_point_equals_lookup() {
        let domain = HashDomain::new(MERKLE_CRH_PERSONALIZATION);

        for _ in 0..TESTRUNS {
            let left = Base::random(thread_rng()).to_le_bits();
            let right = Base::random(thread_rng()).to_le_bits();

            let first = i2lebsp_k(0);

            let input = iter::empty()
                .chain(first.iter().copied())
                .chain(left.iter().by_val().take(L_ORCHARD_MERKLE))
                .chain(right.iter().by_val().take(L_ORCHARD_MERKLE));

            let res1: CtOption<pallas::Point> = domain.hash_to_point_inner(input.clone()).into();
            let res2: CtOption<pallas::Point> =
                domain.hash_to_point_inner_lookup(input.clone()).into();
            assert_eq!(res1.unwrap(), res2.unwrap());
        }
    }

    #[test]
    fn consistent_perm() {
        let domain = HashDomain::new(MERKLE_CRH_PERSONALIZATION);

        for _ in 0..TESTRUNS {
            let left1 = Base::random(thread_rng()).to_le_bits();
            let right1 = Base::random(thread_rng()).to_le_bits();

            let left2 = Base::random(thread_rng()).to_le_bits();
            let right2 = Base::random(thread_rng()).to_le_bits();

            let first = i2lebsp_k(0);

            let input1 = iter::empty()
                .chain(first.iter().copied())
                .chain(left1.iter().by_val().take(L_ORCHARD_MERKLE))
                .chain(right1.iter().by_val().take(L_ORCHARD_MERKLE));

            let input2 = iter::empty()
                .chain(first.iter().copied())
                .chain(left2.iter().by_val().take(L_ORCHARD_MERKLE))
                .chain(right2.iter().by_val().take(L_ORCHARD_MERKLE));

            let perm1 = domain.hash(input1.clone()).unwrap();
            let perm2 = domain.hash(input1.clone()).unwrap();
            let perm3 = domain.hash(input2.clone()).unwrap();
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }
}