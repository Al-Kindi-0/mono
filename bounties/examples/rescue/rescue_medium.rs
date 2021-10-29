use ff::{from_hex, Field};
use zkhash_bounties::{
    fields::{field::Fp, utils},
    rescue_prime::{rescue_prime::RescuePrime, rescue_prime_instances::RESCUE_PRIME_PARAMS_MEDIUM},
};

type Scalar = Fp;

static RANDOM_INPUT: bool = false;

fn main() {
    let params = &RESCUE_PRIME_PARAMS_MEDIUM;
    let rescue = RescuePrime::new(params);

    println!("Rescue Challange medium");
    println!("N = {}", params.get_rounds());

    // insert your solution here:
    let solution: Scalar = from_hex("0x0000000000000000").unwrap();

    let input = if RANDOM_INPUT {
        [utils::random_scalar(true), Scalar::zero()]
    } else {
        [solution, Scalar::zero()]
    };

    let output = rescue.permutation(&input);

    println!("Input  = {:?}", input);
    println!("Output = {:?}", output);

    if output[output.len() - 1] == Scalar::zero() {
        println!("Challenge solved!");
    } else {
        println!("Challenge not solved!");
    }
}
