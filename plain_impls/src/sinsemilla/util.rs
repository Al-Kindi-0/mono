/// Takes in an FnMut closure and returns a constant-length array with elements of
/// type `Output`.
pub fn gen_const_array<Output: Copy + Default, const LEN: usize>(
    mut closure: impl FnMut(usize) -> Output,
) -> [Output; LEN] {
    let mut ret: [Output; LEN] = [Default::default(); LEN];
    for (bit, val) in ret.iter_mut().zip((0..LEN).map(|idx| closure(idx))) {
        *bit = val;
    }
    ret
}
