use super::rescue_prime_params::RescuePrimeParams;

use crate::fields::st::FpST;

use ff::from_hex;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = FpST;

lazy_static! {
    pub static ref MDS3: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0x000000000000000000000000000000000000000000000000000000000000001b").unwrap(),
            from_hex("0x03f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffda").unwrap(),
            from_hex("0x000000000000000000000000000000000000000000000000000000000000000d").unwrap(),
        ],
        vec![
            from_hex("0x000000000000000000000000000000000000000000000000000000000000015f").unwrap(),
            from_hex("0x03f9fffffffffffffffffffffffffffffffffffffffffffffffffffffffffe21").unwrap(),
            from_hex("0x0000000000000000000000000000000000000000000000000000000000000082").unwrap(),
        ],
        vec![
            from_hex("0x0000000000000000000000000000000000000000000000000000000000000db6").unwrap(),
            from_hex("0x03f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffed92").unwrap(),
            from_hex("0x00000000000000000000000000000000000000000000000000000000000004ba").unwrap(),
        ],
    ];
    pub static ref RC3: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0x030b55a63326a9fd7f41acc0d7d0ac16d56b41a2cb6436a1a4208cdccaa91c5d").unwrap(),
            from_hex("0x0015a20e85e4a2522cb9e3c492aa4e0585f5d6e5ab974188067daa47594ccdd3").unwrap(),
            from_hex("0x027753ec9ae7467a393095cbe6662214341d7209f048bf849f3e3bfb247e3f27").unwrap(),
        ],
        vec![
            from_hex("0x0276da3c906ad373bf8f7275dd129d002b39c02f58f6c55e4e064c9453544bc5").unwrap(),
            from_hex("0x01be1d667c315d284729bd0cc4f536f77dc629c592f7634a88a79227025a6638").unwrap(),
            from_hex("0x0276303f2607da2549dbf2e826b8bc15076b1cd22017114526a2f8f2e6d6851a").unwrap(),
        ],
        vec![
            from_hex("0x01fa1071828d84095cdcb78382407e21792f94b8a5f527f6e7c6b8b145a99a6a").unwrap(),
            from_hex("0x0245cfd1470f3d7529363257fb99a185f1ba2a638db0af45eca3bd09d9b83032").unwrap(),
            from_hex("0x02a16181f9fdf734283e742afbac74ff1d92535a88490d103cfc939874d8ea22").unwrap(),
        ],
        vec![
            from_hex("0x015f9db5c45e5a3d9ed2f2e0a073ec7e9b61767890cd0d2b15d2e0fe93ee6097").unwrap(),
            from_hex("0x00ebb6fd31fbae5fce0f7326eb55f3b5858e07ba9c8b1f587a6161ce935ca751").unwrap(),
            from_hex("0x009f0e278eaefd89bf4b05d1e763511ba26bc76fa628fb389e9f830a102867b7").unwrap(),
        ],
        vec![
            from_hex("0x0268e29fb6f2581111ff7e1ba91c0eaf00f334c71957353ec55f155567178f74").unwrap(),
            from_hex("0x00d0ba90e5eaa69b51de0590421a2a0508b822ff1ee560f04fbbd55300c544d2").unwrap(),
            from_hex("0x03e98de3983121fb3bc348d221b8df8d2e49d0663401fb2517f2cb26862ef892").unwrap(),
        ],
        vec![
            from_hex("0x01302414397e1186e6135df77803d74fe69ba908eaf32fbea42b56b42bbe5ae4").unwrap(),
            from_hex("0x0134652331a3b22461d2655c670fac458cbb53fe6d4eef6cf38a8a14bd4ce695").unwrap(),
            from_hex("0x01776b16f9e4603b62cda49bdd2f2bf39873ea62835d43559fd3e50630ebb6a7").unwrap(),
        ],
        vec![
            from_hex("0x012666e58567ead6659491cb6f00dbe2210703ff615e7ec4aeb1c3342ce3ebe7").unwrap(),
            from_hex("0x03ccecfdf4111359b27f30063dc1205a573b26054977a48d80162f3344fae3c7").unwrap(),
            from_hex("0x01e7860342e371165ef75afc73b102d1eca33af8e9f0677835225ca5f2d4a0e6").unwrap(),
        ],
        vec![
            from_hex("0x012a98cfd2fadfa795ef56dc7b5a74700b8fa0abaa082aa59051f112bc842c41").unwrap(),
            from_hex("0x02e8b73c1d2c0fb38e16b6c3441023f136aaf2a89d52ae658cb43360bc90c9b7").unwrap(),
            from_hex("0x022417ee8933d1493085fa862cc46756fb1419760b5a31c87df9911d9b867333").unwrap(),
        ],
        vec![
            from_hex("0x0296639501190ea5e899c5c18d9cac78b43f5032e26c3c9363bc6c64c0467793").unwrap(),
            from_hex("0x0323f4d0a206c6389428d296e0916388783d27c924ee78afba3178a300ab76ea").unwrap(),
            from_hex("0x004fe780bc2eba3222e44e53bec8791fdf4a51d452ebc432a983dcd443bd31d7").unwrap(),
        ],
        vec![
            from_hex("0x02ebd22d13613618e51dd59965bfcfa43155d64f7a6ebe08254564d336e4dd85").unwrap(),
            from_hex("0x00637e60e5795a9d10d4f88711bc24d3ec6392b84ad0e2a490deb0aff1c44fbc").unwrap(),
            from_hex("0x00a2f4fb56f45d9ace0236fb643116d47bbef8996bfe213bd09b6cd3e889bfa4").unwrap(),
        ],
        vec![
            from_hex("0x031d1e52e9e41afa3c1b269f021c47decf8a3dc7b883584fdabb5f02965df171").unwrap(),
            from_hex("0x02ee1bafee2d838462beba2755e15dc4c95536415a317ed6a3a638a554f9944e").unwrap(),
            from_hex("0x01d3b508d3c61e769386c5315fb98d0e33843aaf308a03c1c831b3ab1cd24f11").unwrap(),
        ],
        vec![
            from_hex("0x034490ced4f0daf2f7148e5e1ef21baf097b4dc1bfdfe99f0442787052d6d677").unwrap(),
            from_hex("0x0125d6af33e25a0f8632225c5c7f042ea3f1751da90caa23b1de1f004923088d").unwrap(),
            from_hex("0x00b28a51a0b570e9cd7b8eaf4b33f0ffa568e49d085fdb9a31e156e987271c07").unwrap(),
        ],
        vec![
            from_hex("0x023407a97a29c4cf9ca527b92c86d4c64755fa1713b2c881615bdc043c02aa2f").unwrap(),
            from_hex("0x029971e3a182d3232dca3ea176c7dfd2c148151afbc0c671de89908df3a02bca").unwrap(),
            from_hex("0x02e9a9ee920496c336b0a34851861df506974a79f5dbde0a600f7393d15f6383").unwrap(),
        ],
        vec![
            from_hex("0x024389e7c248d5c178d4630234b5fc8aa933ff8f554d6067b2396e285294eaff").unwrap(),
            from_hex("0x019ebd53ba2cbf9ec9027758a71f390c0b8d4c2fb1a3e834791b153e07ef1c90").unwrap(),
            from_hex("0x03ce59c320593951bf4504244f1eeeecf6624dbf16623accf9447072992dbdb7").unwrap(),
        ],
        vec![
            from_hex("0x0272533f7abaadde46aa3c32f8e42d4acadf488a81538f4700280e9542e67d1c").unwrap(),
            from_hex("0x037aa4f504bd37a6e002f446d2ddefecea012c3fa5888676ed2cab77aaf4c6b7").unwrap(),
            from_hex("0x01e6614a7c1397062f26b8885617cd8efce9f263f276b287a9fb32fb658fc775").unwrap(),
        ],
        vec![
            from_hex("0x03b3f624504df27843631c45f5a1232be9be18783d8306913a31750ee2180767").unwrap(),
            from_hex("0x00beb4427148a3e02fdaae75cb18b7d7cb98a2275b5cabfc0a9e00b91199345b").unwrap(),
            from_hex("0x0059d17d6b6358aaa2f6fed976d52f3c0cd7a83b3266b9a29bc416bdd0c03896").unwrap(),
        ],
        vec![
            from_hex("0x0065b7e5beca7283a3e6bf7e7a1ea9bbce24d215c66f804956fc4ff52b52bd92").unwrap(),
            from_hex("0x025252bc7d5236fa59ffdd4ac0a623a4d90fa6a5cab6cd1e7534de888fd27e6b").unwrap(),
            from_hex("0x0390506b85f297f4640a204ff7fa4f5caf4905f5652d8284c79d3b4fd1147eac").unwrap(),
        ],
        vec![
            from_hex("0x0247d41ba649026749af9d2b5e7752575cad9324f29f025214bcec971eaa1f5e").unwrap(),
            from_hex("0x02a6d67961d2c8140855cdc9f42bd7a7821f870789d04734809e127672c80b77").unwrap(),
            from_hex("0x027ea9dac308a508a36003d4bf77d1c3dfbdf156082bc3b7090a16798db0291b").unwrap(),
        ],
        vec![
            from_hex("0x0156c5a0bcb70e6025544a50a1d226061d06897497b23f72390e9869cbad0f2a").unwrap(),
            from_hex("0x01f64eb8bf610c939a1b0808adee90c5b8aeab5d576aedbb3a6334a70a1c02bc").unwrap(),
            from_hex("0x02be60a37230cf2644107facad363f0b18829c6f9cbf5e52df4ceb1f89193202").unwrap(),
        ],
        vec![
            from_hex("0x03552ee8c1e3441faa995d19d544161b53f8dbbe241bdbbdea4ce3f7b08cd8c6").unwrap(),
            from_hex("0x013bc72dc11b465ccdad150340e05ab4163e5bbcbcfb23abaa0f2c060a3bc118").unwrap(),
            from_hex("0x02d2f5dae8d28a4845ef8148d222e91e99670b083445f98de60f2469693d8026").unwrap(),
        ],
        vec![
            from_hex("0x01a7e556c34ab934f4ab27c46a7fdf184ae382ef6f48543a041cc469c393e5bb").unwrap(),
            from_hex("0x031fc111fc241055f27fd493a16b81e693274aae58e57e63ff6c61e01e4bb6e7").unwrap(),
            from_hex("0x02c6914420fc94b1abed7a75e56cd39faac4e9e5faf8687c641f7a773e843114").unwrap(),
        ],
        vec![
            from_hex("0x002c2762c4a06999b4eb3ce3a7a7c586ccbbeb462fe2882169ef4e4a97a60c6d").unwrap(),
            from_hex("0x01f79cab93cf2475d11a0f5650b8986efef15c349fc8609502921ad415d58304").unwrap(),
            from_hex("0x000df88496b6a38bcffcd18152539750d95ccb6aa7cd2f7f26cc98c13b9e82f1").unwrap(),
        ],
        vec![
            from_hex("0x009ffcfcf75598851997b08a9b03d05d24c759defd8ed43f22ca57cc34c84652").unwrap(),
            from_hex("0x01ec41c978ead905043dedbea56bf5f6ee55fdd9be7f3d7ffec7a573bbb75f96").unwrap(),
            from_hex("0x02f14a23f0dc0acc027463dd5c2e5995dfc8eef5c2d291fce160d6ede6c8c463").unwrap(),
        ],
        vec![
            from_hex("0x016b711a489f7b64a63c158448e1c77eacfd053a7bca9282f997b5c21a417d65").unwrap(),
            from_hex("0x00ffe818f75a487f4d40462a71af82bbbac2355e7131a64546e8c4d8abbe2d98").unwrap(),
            from_hex("0x0211b244341d28a3974a89bdc27063464cc363ea28753e48c6fdf5678d0e27bf").unwrap(),
        ],
        vec![
            from_hex("0x03a9ae5d5103fbb6a4be39f47d72084d795a09aaa180f17ec1a89878c7fe51cb").unwrap(),
            from_hex("0x002803ef7179c0676e69dd056a5f96a64e305c4b0bfa3b2de75f94c2068f1319").unwrap(),
            from_hex("0x019cd5b0dec78d32d1e8e22b2afaca371b164d0eec36fc2e35934df91ce644f7").unwrap(),
        ],
        vec![
            from_hex("0x028e8813a6b991d1bccd37cdc7059d01d2b1ec3fa7e3396388388408dd57247d").unwrap(),
            from_hex("0x023cb2848f9b84d5a551030a32ca9bd1c9b0e0af09acc9bdf00f214583e32d57").unwrap(),
            from_hex("0x019f84d3b8ca2f5726cb28be191b7d1ca6c62b46c5454089226a6005d0d414ca").unwrap(),
        ],
        vec![
            from_hex("0x005a8f7c640dcd59e7ffb6ebca4b9a59e83b033d052527a030d4c5ab3933c0c1").unwrap(),
            from_hex("0x00c896a41e3e4611ee312fe9750727181d00a6b13234a9256eaf4b8d9f8762f0").unwrap(),
            from_hex("0x02ade9f09680bfbee2bec3cce9bdf67bf328102ce963dd6a3667a2dc06b32af9").unwrap(),
        ],
        vec![
            from_hex("0x00eaf54e2930c6b0362186b43b91f7ef00f994b38d3d59cdbaddf47dc7528fd0").unwrap(),
            from_hex("0x0321f660e39fb04c50cb9612c0a14576fc60dada1abfa6eab479b4fefe3cea3c").unwrap(),
            from_hex("0x03957c77885fe30a15ba2a30cf0a6073164ff8f46de4181ab0a38999b2683bfe").unwrap(),
        ],
        vec![
            from_hex("0x03a205569f0ffb955744635aafbab2adccce0fb573f7bd08d1b2fbe76aa7e625").unwrap(),
            from_hex("0x0240b2d450acd64195ba91fae37691cf307f57d9a6247c63e3bfc7341ab9e5ea").unwrap(),
            from_hex("0x004b2ac1983dfcedff20797065421df1489803334fede41a77720fd59c21e4bf").unwrap(),
        ],
        vec![
            from_hex("0x03cefe5ccddbdc2e5b4ae27eb4f51259f36fe2a7db8c064e622b946a9a48d0d0").unwrap(),
            from_hex("0x00e1b2f1ac1b12f2980e824e95ec42a15d27ac8c9521eac62cd917ad89cb944f").unwrap(),
            from_hex("0x0171adc778a0f128a97ef57e1f17153b47fec7b4d422c88452fabded85f46905").unwrap(),
        ],
        vec![
            from_hex("0x032aea9cb8c83eef31eb872bf510cd48bb10c747b92dafee96b13735e6329783").unwrap(),
            from_hex("0x03003a24cfcc8ce303dc123eeeeaa18773d3ce6e82d8dc0dd4863a4e8d35fc96").unwrap(),
            from_hex("0x014485c92a86453735964299412deb54d5f8b619693ecce1d10f10b187561301").unwrap(),
        ],
        vec![
            from_hex("0x00525659ac4e172aa8e605fe453eadddf3946cfdc4c23308a46d7f9d6e703550").unwrap(),
            from_hex("0x01a5925666e33813d618494f0ad4d791bd7598f15c5cc658f4a377e265764e8c").unwrap(),
            from_hex("0x018330e714ca4f02ea1a4f41cdfc8816127399ecf161926e8148e7c2e17eec62").unwrap(),
        ],
        vec![
            from_hex("0x00648c61c689136e4fe4b4fa596877466282d19f2f2f44e267192a0e652f7ccc").unwrap(),
            from_hex("0x00f1f9ec81c982c0de6382f9778a45c7b58aad1e79337f93c1132ee39c5a2287").unwrap(),
            from_hex("0x00abc5a11c38bbac867754138a4e6f85ab643f664c818847ed341bc702c6126b").unwrap(),
        ],
        vec![
            from_hex("0x00aa4c2f0d9d3621d1a5b317770b700b9c2df01f13c248467ea312f4d31dc7c1").unwrap(),
            from_hex("0x02d11bed5b6f6ed20ef5fbd6129425a5bfb3949c3563209d6e127c553e52cc9f").unwrap(),
            from_hex("0x029e4ff0f088405aba57f15f01d17e26da64c1cdd22d63e23390f7d5d2d52de5").unwrap(),
        ],
        vec![
            from_hex("0x01d1cd9cdd619ec51754b15b081c5cc1f06f8ba916ae514722bce154607c3c81").unwrap(),
            from_hex("0x036acc75ede39d8ce78e7abe45c03001f759b3310d776378b0c81c221791b8b0").unwrap(),
            from_hex("0x00da7c1c15d42c1479da49a35f81e483d97fdfd982d348eeba4ed443f808a9fa").unwrap(),
        ],
        vec![
            from_hex("0x02fc7050bd95c7a878fbbda02f7a9fbf366330e5714c582f21a290d6ccbced2b").unwrap(),
            from_hex("0x03d927539a97a3d925241fed2fc9dea9cc57379b2a5569bbde1ff738a69bab8f").unwrap(),
            from_hex("0x0179342934e52633fdde3e2f3d7c3235e1dc0e8f8250bfc3b67e7d07c53b2626").unwrap(),
        ],
    ];
    pub static ref RESCUE_PRIME_ST_PARAMS: Arc<RescuePrimeParams<Scalar>> =
        Arc::new(RescuePrimeParams::new(
            3,
            3,
            [
                0xaaaaaaaaaaaaaaab,
                0xaaaaaaaaaaaaaaaa,
                0xaaaaaaaaaaaaaaaa,
                0x2a6aaaaaaaaaaaa,
            ],
            18,
            &MDS3,
            &RC3
        ));
}