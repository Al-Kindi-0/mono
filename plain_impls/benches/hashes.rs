use std::iter;

use blake2::{Blake2b, Blake2s};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use group::ff::Field;
use group::ff::PrimeFieldBits;
use pasta_curves::pallas::Base;
use random::thread_rng;
use random::Rng;
use sha2::{Digest, Sha256};
use sha3::Sha3_256;
use zkhash::pedersen_hash::pedersen_hash::pedersen_hash;
use zkhash::pedersen_hash::pedersen_hash::Personalization;
use zkhash::sinsemilla::sinsemilla::{
    i2lebsp_k, HashDomain, L_ORCHARD_MERKLE, MERKLE_CRH_PERSONALIZATION,
};

fn sha256(c: &mut Criterion) {
    let input = b"hello_world";

    c.bench_function("SHA256 Hash", move |bench| {
        bench.iter(|| {
            let hash = Sha256::digest(black_box(input));
            black_box(hash)
        });
    });
}

fn sha3_256(c: &mut Criterion) {
    let input = b"hello_world";

    c.bench_function("SHA3-256 Hash", move |bench| {
        bench.iter(|| {
            let hash = Sha3_256::digest(black_box(input));
            black_box(hash)
        });
    });
}

fn blake2s(c: &mut Criterion) {
    let input = b"hello_world";

    c.bench_function("Blake2s Hash", move |bench| {
        bench.iter(|| {
            let hash = Blake2s::digest(black_box(input));
            black_box(hash)
        });
    });
}

fn blake2b(c: &mut Criterion) {
    let input = b"hello_world";

    c.bench_function("Blake2b Hash", move |bench| {
        bench.iter(|| {
            let hash = Blake2b::digest(black_box(input));
            black_box(hash)
        });
    });
}

fn sinsemilla(c: &mut Criterion) {
    let domain = HashDomain::new(MERKLE_CRH_PERSONALIZATION);

    let left = Base::random(thread_rng()).to_le_bits();
    let right = Base::random(thread_rng()).to_le_bits();

    let first = i2lebsp_k(0);

    let input = iter::empty()
        .chain(first.iter().copied())
        .chain(left.iter().by_val().take(L_ORCHARD_MERKLE))
        .chain(right.iter().by_val().take(L_ORCHARD_MERKLE));

    c.bench_function("Sinsemilla Hash", move |bench| {
        bench.iter(|| {
            let hash = domain.hash(black_box(input.clone()));
            black_box(hash)
        });
    });
}

fn pedersen(c: &mut Criterion) {
    let mut rng = thread_rng();
    let personalization = Personalization::MerkleTree(0);
    let input: Vec<bool> = (0..510).map(|_| rng.gen::<bool>()).collect();

    c.bench_function("Pedersen Hash", move |bench| {
        bench.iter(|| {
            let hash = pedersen_hash(black_box(personalization), black_box(input.clone()));
            black_box(hash)
        });
    });
}

fn criterion_benchmark_hashes(c: &mut Criterion) {
    sha256(c);
    sha3_256(c);
    blake2s(c);
    blake2b(c);
    sinsemilla(c);
    pedersen(c);
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark_hashes
);
criterion_main!(benches);
