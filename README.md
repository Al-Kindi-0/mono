# Plain Performance Comparison of different Hash Functions for ZKP

This repository contains Rust implementations of different hash functions for Zero-Knowledge applications.

## Hash Functions

The following hash functions are already implemented:

- [ReinforcedConcrete](https://todo)
- [Poseidon](https://eprint.iacr.org/2019/458.pdf)
- [Rescue](https://eprint.iacr.org/2019/426.pdf)
- [Rescue-Prime](https://www.esat.kuleuven.be/cosic/publications/article-3259.pdf)
- [Feistel-MiMC](https://eprint.iacr.org/2016/492.pdf)
- [Pedersen-Hash](https://zips.z.cash/protocol/protocol.pdf#concretepedersenhash), code extracted from [Zcash](https://github.com/zcash/librustzcash)
- [Sinsemilla](https://zips.z.cash/protocol/protocol.pdf#concretesinsemillahash), code extracted from [Orchard](https://github.com/zcash/orchard)

We also benchmark against various classical hash algorithms implemented in [RustCrypto](https://github.com/RustCrypto/hashes).

We instantiate the finite-field permutations (ReinforcedConcrete, Poseidon, Rescue, Rescue-Prime) with a statesize of three field elements in a sponge with one field element reserved as the capacity. Feistel-MiMC always has a statesize of two, which is why one can only absorb one field element per permutation call when instantiated in a sponge.

## Benchmarks

Here we give benchmarks for hashing input sizes of 512-bit (i.e., two field elements for the used prime fields). We,  thus,  benchmark  one  permutation call for all symmetric hash functions, except for Feistel-MiMC for which we require two. All benchmarks where obtained  on  a  Linux  Desktop  PC  with  an  Intel  i7-4790 CPU (3.9 GHz) and 16 GB RAM using stable Rust version 1.53 and the `target-cpu=native` flag. Time in ns.

| Hash               |         |  BN     | BLS     | ST      |
|--------------------|--------:|--------:|--------:|--------:|
| ReinforcedConcrete | -       |   3 284 |   3 265 |   1 032 |
| Poseidon           | -       |  17 598 |  18 174 |  17 320 |
| Rescue             | -       | 415 230 | 446 980 | 359 510 |
| Rescue-Prime       | -       | 362 870 | 391 560 | 294 660 |
| Feistel-MiMC       | -       |  33 800 |  35 847 |  28 594 |
| Sinsemilla         | 131 460 | -       | -       | -       |
| Pedersen-Hash      |  39 807 | -       | -       | -       |
| SHA-256            |   366.5 | -       | -       | -       |
| Blake2b            |   245.1 | -       | -       | -       |
| Blake2s            |   219.5 | -       | -       | -       |
| SHA3-256           |   392.3 | -       | -       | -       |
