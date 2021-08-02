use ff::PrimeField;

#[derive(Clone, Debug)]
pub struct PoseidonParams<S: PrimeField> {
    pub(crate) t: usize, // statesize
    pub(crate) d: usize, // sbox degree
    pub(crate) rounds_f_beginning: usize,
    pub(crate) rounds_p: usize,
    pub(crate) rounds_f_end: usize,
    pub(crate) rounds: usize,
    pub(crate) mds: Vec<Vec<S>>,
    pub(crate) round_constants: Vec<Vec<S>>,
    pub(crate) opt_round_constants: Vec<Vec<S>>, // optimized
    pub(crate) w_hat: Vec<Vec<S>>,               // optimized
    pub(crate) v: Vec<Vec<S>>,                   // optimized
    pub(crate) m_i: Vec<Vec<S>>,                 // optimized
}

impl<S: PrimeField> PoseidonParams<S> {
    pub fn new(
        t: usize,
        d: usize,
        rounds_f: usize,
        rounds_p: usize,
        mds: &[Vec<S>],
        round_constants: &[Vec<S>],
        opt_round_constants: &[Vec<S>],
        w_hat: &[Vec<S>],
        v: &[Vec<S>],
        m_i: &[Vec<S>],
    ) -> Self {
        assert!(d == 3 || d == 5);
        assert_eq!(mds.len(), t);
        assert_eq!(rounds_f % 2, 0);
        let r = rounds_f / 2;
        let rounds = rounds_f + rounds_p;

        PoseidonParams {
            t,
            d,
            rounds_f_beginning: r,
            rounds_p,
            rounds_f_end: r,
            rounds,
            mds: mds.to_owned(),
            round_constants: round_constants.to_owned(),
            opt_round_constants: opt_round_constants.to_owned(),
            w_hat: w_hat.to_owned(),
            v: v.to_owned(),
            m_i: m_i.to_owned(),
        }
    }
}
