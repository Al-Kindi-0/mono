use ff::{from_hex, Field};
use zkhash_bounties::{
    feistel_mimc::{feistel_mimc::FeistelMimc, feistel_mimc_instances::FM_PARAMS_MEDIUM},
    fields::{field::Fp, utils},
};

type Scalar = Fp;

static RANDOM_INPUT: bool = false;

fn main() {
    let params = &FM_PARAMS_MEDIUM;
    let feistel_mimc = FeistelMimc::new(params);

    println!("FeistelMimc Challange medium");
    println!("r = {}", params.get_rounds());

    // insert your solution here:
    let solution: Scalar = from_hex("0x0000000000000000").unwrap();

    let input = if RANDOM_INPUT {
        [utils::random_scalar(true), Scalar::zero()]
    } else {
        [solution, Scalar::zero()]
    };

    let output = feistel_mimc.permutation(&input);

    println!("Input  = {:?}", input);
    println!("Output = {:?}", output);

    if output[output.len() - 1] == Scalar::zero() {
        println!("Challenge solved!");
    } else {
        println!("Challenge not solved!");
    }
}
