use tplonk::{
    fields::{bls12::FpBLS12, utils},
    reinforced_concrete::{
        reinforced_concrete::ReinforcedConcrete, reinforced_concrete_instances::RC_BLS_PARAMS,
    },
};
type Scalar = FpBLS12;

fn perm_example() -> [Scalar; 3] {
    let rc = ReinforcedConcrete::new(&RC_BLS_PARAMS);

    println!("Sample input");
    let input: [Scalar; 3] = [
        utils::random_scalar(true),
        utils::random_scalar(true),
        utils::random_scalar(true),
    ];

    println!("Permutation");
    rc.permutation(&input)
}

fn hash_example() -> Scalar {
    let rc = ReinforcedConcrete::new(&RC_BLS_PARAMS);

    println!("Sample input");
    let input1 = utils::random_scalar(true);
    let input2 = utils::random_scalar(true);

    println!("Hash");
    rc.hash(&input1, &input2)
}

fn main() {
    println!("ReinforcedConcrete BLS12 Example");
    perm_example();
    hash_example();
}
