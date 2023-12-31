use super::rescue_prime_params::RescuePrimeParams;

use crate::fields::bn256::FpBN256;

use ff::from_hex;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = FpBN256;

lazy_static! {
    pub static ref MDS3: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0x000000000000000000000000000000000000000000000000000000000000007d").unwrap(),
            from_hex("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593efffff66").unwrap(),
            from_hex("0x000000000000000000000000000000000000000000000000000000000000001f").unwrap(),
        ],
        vec![
            from_hex("0x0000000000000000000000000000000000000000000000000000000000000f23").unwrap(),
            from_hex("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593efffedb9").unwrap(),
            from_hex("0x0000000000000000000000000000000000000000000000000000000000000326").unwrap(),
        ],
        vec![
            from_hex("0x000000000000000000000000000000000000000000000000000000000001898e").unwrap(),
            from_hex("0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593effe2722").unwrap(),
            from_hex("0x0000000000000000000000000000000000000000000000000000000000004f52").unwrap(),
        ],
    ];
    pub static ref RC3: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0x241214b64e37a42dddc49216b6433fe75e4af3533a8c8961def18b459420ce96").unwrap(),
            from_hex("0x149e9522e80164b39561a6d532ed480ddb16db399fce8f2b72c8640bed14edd8").unwrap(),
            from_hex("0x16e6151eb7f6065df49647b709fcde486776be1e372155e42ce9c91b49342af3").unwrap(),
        ],
        vec![
            from_hex("0x0b29463b35fc98ca03baae98f5d4f251d38e091fa179fbe1f10e77e0f46399cd").unwrap(),
            from_hex("0x1a892f66364b75798cebe8e3ef3bf830f85ecb833f3b1023f4e90d2fc67a88d9").unwrap(),
            from_hex("0x27af8ef16eb7a0535a73aaa4273ea0811b95b6f288e3b18e91ea29857a35f4df").unwrap(),
        ],
        vec![
            from_hex("0x03546ef6134d6bcfac31bfdcc211836203a559b0e04b314ff40642b72b1f22b9").unwrap(),
            from_hex("0x15392eceb2d870dcedce619bd4a4baf5140dfe1390fcb22d89787f5631f4756d").unwrap(),
            from_hex("0x238fa99e483edf2d219d37dcf824e86ed4fc7b384a4fecc8a0f5c8109fa0d0e7").unwrap(),
        ],
        vec![
            from_hex("0x00bf1355bb7cfb01c74f7188922e6120a8c24b612471653b59ac1aff07d44f46").unwrap(),
            from_hex("0x16af47292baf23e76016f26496f6c73d3c38f4b1791f3b8762f7f15b89acc9a3").unwrap(),
            from_hex("0x105a902eafac24043e91c89164e510e1d7d1948c5660b56c7f6f7672dbe60b75").unwrap(),
        ],
        vec![
            from_hex("0x0d7bc0cc3063a9d7c2b85c953ac460a79e70b3ac3fef6f7952658a59cce56cce").unwrap(),
            from_hex("0x205aa50ee2dc005f22a93fd5070636e20406843871f8eb4bfce2845e060df6c1").unwrap(),
            from_hex("0x0115db7f2494ba498f5168807e455811e55f92d6b48117dd75791b94b0d2e42f").unwrap(),
        ],
        vec![
            from_hex("0x116ee19eba3b6f6f24d41133d4cc4d6b22d940b45256d7fb55f7a2bbc56fe585").unwrap(),
            from_hex("0x1def0328023519e98741d1bf42429a1ddf9853579b3a59424184c12d0137c33c").unwrap(),
            from_hex("0x005746d2203f013e44ac7cd2ca4f7025881a574932e81ee15125558032ca2d9b").unwrap(),
        ],
        vec![
            from_hex("0x2c7c456a9d460f23f299aa325599598947fb9289a8d9743efa621cea723cfb8b").unwrap(),
            from_hex("0x04f971670f113c12f3e59a3bf8c32a9d71526bd9ff1ca7933bcf42f57cbe9142").unwrap(),
            from_hex("0x1ba8e88c0c59e257fce11428c76c7d5d35b5a773d6c7da2609869306e43d0959").unwrap(),
        ],
        vec![
            from_hex("0x00762d2d87b03a76a8851f6f5a69d3a6078a8be5028f962098da7be06012a7bc").unwrap(),
            from_hex("0x23d7c7a4017398ef348dbb6b4d9f531b0757ad6050704bc7641be4492cf3ca0f").unwrap(),
            from_hex("0x0a9cccd695ee8ad147aa245d3d7a30cfc0e0ab6072910d02cc00cf978ca0a89e").unwrap(),
        ],
        vec![
            from_hex("0x0a0f24024b1fb9b5afd2997d5e3f1793f99e3f0bccb54dadefddd77b7a42a058").unwrap(),
            from_hex("0x024e7d70b40332e5e0d5c790f244ce1685c3062776a792d6144385d1031c6a76").unwrap(),
            from_hex("0x0f3e9d716963356b6f4d6c59b5ba5ddccbe583b556a5466e6a8f2d444191542b").unwrap(),
        ],
        vec![
            from_hex("0x1fb906fad59abc852df6bad6e47237de825ec36ba13efc46f29a4c8dd680bcaa").unwrap(),
            from_hex("0x282fe85ec5d4b5bf5ac6dce1e237e078d107b69b35fe263e295b41a707cb9b42").unwrap(),
            from_hex("0x0dee9f78d30ebbcabc8dc0ef1bda1e0a665c218fc9e1d6a1446cc9c8a765015a").unwrap(),
        ],
        vec![
            from_hex("0x18d28d7bea35db5e6ddacd576186bc4aae4485f623b49fe505a94054ba75ef0e").unwrap(),
            from_hex("0x0cc59f4f8d39b3ac4f2567bc56eb2c7e8da5550a0c04818b42fe695f8f285c20").unwrap(),
            from_hex("0x1f1ef239cea48c9aaafcb216b0e08e5fd68cbea8eda24235a6f0ce2c85609659").unwrap(),
        ],
        vec![
            from_hex("0x19832475410c38053d6b7085501edf215741207e4bc9548afe8b4b179b9fa253").unwrap(),
            from_hex("0x28f1800567daedaa3673eaa304d90334d616dcc9b6a093df07fab20251d2f27a").unwrap(),
            from_hex("0x293042c65e37a4efb3190692bc75c5470076513c77b87b9d3535c10f1c5ed68f").unwrap(),
        ],
        vec![
            from_hex("0x16019c9451b62d42177d2cacd260a15f0de9cdc9ccb26a892bb8d37ac61ce9bf").unwrap(),
            from_hex("0x2323a90bb17a61acbca2205486b44b706cc90fcb4e9900d2970f0df02575c553").unwrap(),
            from_hex("0x000c38d85cf32503c63b8ac156492b25f550ad2afddfb92c9c31ac4b6603e304").unwrap(),
        ],
        vec![
            from_hex("0x2c69e902753e9b71445f40582287929e6379737b578a73b8b6af949d775d880e").unwrap(),
            from_hex("0x046305445d6def7ea13e73f364b5bb76a2480fcaf9e806bd21dc414ec11b7e48").unwrap(),
            from_hex("0x189b2620678f5309ee12fff3424a7d65b75a2f674d53da1d594266f477afe57f").unwrap(),
        ],
        vec![
            from_hex("0x1c9cb3cad66a96d4f9131760344e7093ae4358bac852f6e352ee345e0dcdf684").unwrap(),
            from_hex("0x11cb60a422f7ffcf0a8de27dddd490b6fd93606c37dffc6e8aea256c157fae69").unwrap(),
            from_hex("0x23ff8be08521aaa5a6ea6c7fd2c7526afca282b354ee7a559099190e072e3ce4").unwrap(),
        ],
        vec![
            from_hex("0x0d4ade548e38a7c4a1976be0cf50cac82e37f202e99413930834a5e117b34276").unwrap(),
            from_hex("0x14bc69cef73fe0bb617b6d21dd01cd7f635b169a8b975b8562d9d8460c1aade8").unwrap(),
            from_hex("0x0db842e9b71b286915efa0de5e03f8a0378b72b7a71f0c2135e79866d3d6f528").unwrap(),
        ],
        vec![
            from_hex("0x1ab35e2f964a0c9641ae01e04747e2a686c76da44ce42579c58c38235ad2eb0c").unwrap(),
            from_hex("0x18de353617b2891c392d9f3b6386d74a81f5c4468eebfe8c73200114972fe5b1").unwrap(),
            from_hex("0x180b471ce6b043a9401cdd596456af7c67b2b800474bf4fc6932e2edbc62cfbf").unwrap(),
        ],
        vec![
            from_hex("0x23153ccd41fb458e2f33c20b8ce49b8985836a26bb06c39f0a96d6bbbc0301ba").unwrap(),
            from_hex("0x1416013abc7d9b53aef83185611f5617aba83bfdda11dd489d77cd2012e8a8a7").unwrap(),
            from_hex("0x220a789dc01b985c3a137384c37d0b5ad86b7f07f6e0224c66fc0798e9f6459c").unwrap(),
        ],
        vec![
            from_hex("0x134dc23093822f920d9c9301b363b224d4fe4f977e11e7a1393244cdfc88ce1a").unwrap(),
            from_hex("0x0a2b2d1d9cc5de93ea90b33dcdbf156c000b753b60db180823e717d2c49d6910").unwrap(),
            from_hex("0x0091ceffd5b51b15c3b608dd743c9a36eb2ffb6ef374011eb5d0ed60f1ff2b49").unwrap(),
        ],
        vec![
            from_hex("0x024de554062063c0168c82ddb650ea09b415c7405d5f3b1430f6eeec0ad3fb0d").unwrap(),
            from_hex("0x2c450f23635a10c72b8fc7f36750643a42f62453cd501c83c6c90d16d7eefc57").unwrap(),
            from_hex("0x2f14c4092eb0a874c85ae64b4d18bba47970ca3da4d422629ceb040e62b14096").unwrap(),
        ],
        vec![
            from_hex("0x2e2561cf8692bdcb2136b5038feca8b05e379a4c6ad4d0f5b4b5af8dc94dc1f7").unwrap(),
            from_hex("0x17310c87b9b20d078bb4ea19756cd049afb5dec9734e9745dccd521b007258e1").unwrap(),
            from_hex("0x0093cb39757463eb403a16afbec58d3fe5bb0db9a0e69b39e23272fe5420b828").unwrap(),
        ],
        vec![
            from_hex("0x0b057a4cb37d03a96cdb20c1f9a96eea600fcd17d2479b228fe9e7ea4befd3c6").unwrap(),
            from_hex("0x189552e5eb3ac601a687cae3675cd9c2f72b1a41bcaee2835c03152078590107").unwrap(),
            from_hex("0x0518b70350baac601b679aca4238937870f9578567ce3696f81bee468fd7024e").unwrap(),
        ],
        vec![
            from_hex("0x0bb790131ed126376809f10b87f7efcd21502e65311c2960f54b0a6953446419").unwrap(),
            from_hex("0x01957d4149e870b8124d9da7c079dbe78471228feba9d21b755d2b76b74e2d0d").unwrap(),
            from_hex("0x049d840bbb1007263ec4103d9c8fefe67581bfd5cedf6a85d2c1b13613d99b87").unwrap(),
        ],
        vec![
            from_hex("0x16663bd42a4d96e3b69edcf1a11950b22cee06402894fdf330128952e31dd397").unwrap(),
            from_hex("0x25299e3a923fc0c38ec4d2421077c7e547d9cbd7f9fc89d2915de94b188417c4").unwrap(),
            from_hex("0x2a238002a34a8c72b397392f399af21e0e0f7fcc05506e6b874dd3f70bba4b3b").unwrap(),
        ],
        vec![
            from_hex("0x0dc6f6cc8d865f25ed467bbf95bbd398dff0d10ae52e2a30be6c8e82344b3802").unwrap(),
            from_hex("0x2a72d90ddf392777d2a1977eeb51cb50aee2e9e8b7b7d4c94fb141de94b1ed69").unwrap(),
            from_hex("0x00160b8013f8d967f070ff7f978763d30f9f50208558cec97a5e3bb4af933a2d").unwrap(),
        ],
        vec![
            from_hex("0x12fa0490ff006e46e16e8aaccba07cd2ba266bfff29bc68449e007349a888207").unwrap(),
            from_hex("0x23ddf606bb111b9b21afc99904fa6f75fe37fce69a7e1f7fd5e3d1292ba7e7e0").unwrap(),
            from_hex("0x29b2ec689b6f2ed0dbf269c777609432b038a431d432de96688047edc088c1a4").unwrap(),
        ],
        vec![
            from_hex("0x2061ba6a4ad4076d895e99f6210c642f745e1e0c103d5407dab83cda773b7fc3").unwrap(),
            from_hex("0x00e536b883f7c592c1f6d648bdffa05fd559f399a0a415c7620c16ee8583bbd1").unwrap(),
            from_hex("0x2d782ff8b4ff168929034808ce0a8b6493d302444e7777d4aff02f55a3b56769").unwrap(),
        ],
        vec![
            from_hex("0x2f6ead5d361bfd4e9986132970afcb85b2f4b4538a6f70bbc75e0f0034c8761a").unwrap(),
            from_hex("0x2ac8b859e3deff5e0036d4bb0f0393bc9311336ffea5122cf5e380ebb82d3b55").unwrap(),
            from_hex("0x0a1e0608a08cdf1baa58db694c6f73b8d6d598ea447dd6bded6fea2470b38d0a").unwrap(),
        ],
    ];
    pub static ref RESCUE_PRIME_BN_PARAMS: Arc<RescuePrimeParams<Scalar>> =
        Arc::new(RescuePrimeParams::new(
            3,
            5,
            [
                0xcfe7f7a98ccccccd,
                0x535cb9d394945a0d,
                0x93736af8679aad17,
                0x26b6a528b427b354,
            ],
            14,
            &MDS3,
            &RC3
        ));
}
