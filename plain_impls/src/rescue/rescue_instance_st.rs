use super::rescue_params::RescueParams;

use crate::fields::st::FpST;

use ff::from_hex;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = FpST;

lazy_static! {
    pub static ref MDS3: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0x000000000000000000000000000000000000000000000000000000000000001b").unwrap(),
            from_hex("0x000000000000000000000000000000000000000000000000000000000000015f").unwrap(),
            from_hex("0x0000000000000000000000000000000000000000000000000000000000000db6").unwrap(),
        ],
        vec![
            from_hex("0x03f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffda").unwrap(),
            from_hex("0x03f9fffffffffffffffffffffffffffffffffffffffffffffffffffffffffe21").unwrap(),
            from_hex("0x03f9ffffffffffffffffffffffffffffffffffffffffffffffffffffffffed92").unwrap(),
        ],
        vec![
            from_hex("0x000000000000000000000000000000000000000000000000000000000000000d").unwrap(),
            from_hex("0x0000000000000000000000000000000000000000000000000000000000000082").unwrap(),
            from_hex("0x00000000000000000000000000000000000000000000000000000000000004ba").unwrap(),
        ],
    ];
    pub static ref RC3: Vec<Vec<Scalar>> = vec![
        vec![
            from_hex("0x02a44b57f6e9b0e2e8b817dac3698c80fb839be67095dff247fa5c7e9dbcedcf").unwrap(),
            from_hex("0x011a45e31139695c94ba33cd73c90040fd0c31b012d02fd1501a72b3f9001ade").unwrap(),
            from_hex("0x02486100be21d1eba5a331ff3b7a3a85698329ed00647810bf7c9ce76060c905").unwrap(),
        ],
        vec![
            from_hex("0x03813ecf7314b02d0fdbf321538d9f1c63585ad5cd0ff65b682a2614b82f9a7a").unwrap(),
            from_hex("0x0148bfee115b2d6d2d5f29a108dbc040dec9a09b81ab754626a3156b1cfe6fb4").unwrap(),
            from_hex("0x00162fc70527380d716ecd150155e97f6b0d98ac7c09a6fa2a41de0577e47f55").unwrap(),
        ],
        vec![
            from_hex("0x02701b5ba8ae980a5487fb12c98085734a7aec600b937b9f135685ad47be680e").unwrap(),
            from_hex("0x01345f15c574482f408946ecd70f6dfc0aace571fe0fd6db7168e59d42dc4d5f").unwrap(),
            from_hex("0x00cabfe5415f9b71c39dcfee3474e4fb6567b039740723796b96d2d9034bda84").unwrap(),
        ],
        vec![
            from_hex("0x027ddea62a0f0bfc487016025effa10f023310d441b034ed029cf25c2b8ae0f7").unwrap(),
            from_hex("0x0311eacb597d16d0f69553d83d5ad717df3c8b770c99cd129b7d7db31c1ef2b9").unwrap(),
            from_hex("0x00f066e86374a43b8a4d7bb1de4ada40357236461c76a15dfd3244acae901c27").unwrap(),
        ],
        vec![
            from_hex("0x0326279d2665b313303b0f491176287179742ce4a3199549c83cb582a8716efc").unwrap(),
            from_hex("0x00eedc5b82f3f5f5bdf2f79e2f67fb02fb23dae43277bd81d71ed4e44637cc2d").unwrap(),
            from_hex("0x0272938073f9eff03592a1f421562c16bc08e1c0f003b11a4c0d52aecf4de1ef").unwrap(),
        ],
        vec![
            from_hex("0x00e601e794bb231f77bc40ea237d4810524536921cff550168cc61ea01d01434").unwrap(),
            from_hex("0x03bd1dea4ac1e9f6d426700c25a95eb87f8c36d72b29e0db8ac409a46508ffdf").unwrap(),
            from_hex("0x007760cab34b042b78679ec1d298fc04e07ae07c6acd64c063ea462fa3f578ff").unwrap(),
        ],
        vec![
            from_hex("0x01e00b4814cf7e4d86ef13a585feced5a6ece7a7eb92da3ccb8090a07dafd2d7").unwrap(),
            from_hex("0x01002c9abfbe6f9438a61d37a9755728ed3ff688e93bcae0a252c42b293b20ea").unwrap(),
            from_hex("0x0360b609c0699270d44c4563a309276f0b7d266251cb7f1e041971abce6be672").unwrap(),
        ],
        vec![
            from_hex("0x02d12939e3bf4314aff9eb5ce05e9445541e4f4e67cb4520ae29e93f7f7a5a95").unwrap(),
            from_hex("0x00e620cc01cc196d34266bbb4e7ab49d5618bcac8ca9ac496dba1bc502409c8b").unwrap(),
            from_hex("0x018108442bbc03ab887cf2730e51cd29673fa329f75871c792edccd0b38a4233").unwrap(),
        ],
        vec![
            from_hex("0x016488e512fc5a8c5217e79dd2e01e8087085302f1b70b65e3776d0d41f54622").unwrap(),
            from_hex("0x007331f97cf5e9d54e9f0a200216c090e5f885238a75ca56568042b6e4edae1d").unwrap(),
            from_hex("0x010448e1b85e40f2e4f315c061413e9f81c9226c8d2b60c4c9cea113ea89a4be").unwrap(),
        ],
        vec![
            from_hex("0x0135990cfcbf9e4b1e4053021129f824dad9530d1d1ece24f69460ada0b2cb65").unwrap(),
            from_hex("0x02a0d7730ca483aa171e9dfd3937c464e75e435fd585532d525ecb2a389b4e22").unwrap(),
            from_hex("0x036be72a677fdd658864451bd3f17b1eb37eba6618b637d52f54f2aabf03aeae").unwrap(),
        ],
        vec![
            from_hex("0x014df48f4787305ce815dd6121c8e3cd741aa13237dac0912770bf5248f85f9a").unwrap(),
            from_hex("0x01a817d243c27cb28a80f8bbfa5210b3f9bc074adadc572a1da934ca2f014c83").unwrap(),
            from_hex("0x009c3ee108a5cc2771e1c83a18dfa21b2fa242967b5ce09ca2d846ff1cf7c4d0").unwrap(),
        ],
        vec![
            from_hex("0x00e2d5d834c2074d50bd0b952e7ede5c2d38d3eaa43403ccb10b3754e2c6d093").unwrap(),
            from_hex("0x00d6bac7994a06ea6ec0603e51669e7e8c22c35d46189a16e9f9b4da06d3b0b6").unwrap(),
            from_hex("0x00cf9e97c3099e689ba3edc8347cf637ebe96ab1099f7811aa358f30ce3bac20").unwrap(),
        ],
        vec![
            from_hex("0x0331a6f4962f71450f7c7db4efeb5a28ea737c6a0894e279a66bb29cef5adfe1").unwrap(),
            from_hex("0x008ea33fc1436f4df4281bc19e697f55e218d2619958448adbaf480465a2ab4e").unwrap(),
            from_hex("0x038d56eb51577107ef576c9d5fbac44b811d7c733bb849b6926acd1c40f35cd2").unwrap(),
        ],
        vec![
            from_hex("0x011e02ebf128d1eaa566a29673c5001d5458271a428eec1fa635fdcaaa976e52").unwrap(),
            from_hex("0x027cb47132310ade4eba490fed90cd67e08b089a5fe18dba826e0dfa9a744072").unwrap(),
            from_hex("0x0387fcd6556435f52a82b00f94c995cbeeb16290c466e40b125eeecab4fc71cb").unwrap(),
        ],
        vec![
            from_hex("0x03c27bc38e356c1be1f19423c5e7074b790f269c961a0d8c9938e40d27596420").unwrap(),
            from_hex("0x016e87f2f7fe456067141d0c5288425fd4a9c1b82ab82cf3e0cf098ae1500a16").unwrap(),
            from_hex("0x03aaa53a0f825486551d0f5bb0a3eebf1429de2f6299a73ddf5b4172bb3f7822").unwrap(),
        ],
        vec![
            from_hex("0x02db95a729b2247be47989aa6f7ca7d2214326fcd3d0405d2e58b6458f6e3197").unwrap(),
            from_hex("0x034fd13929dc1867f1a7d3116ba8eb8afd488e8e8a0cadeca57536def02a37a1").unwrap(),
            from_hex("0x00caeda58a453bf7c06e21a7f7d9db18fddf6b7eb453f6fbef568dacba0ed00b").unwrap(),
        ],
        vec![
            from_hex("0x0128cfd9e3c5c503db72c5bddc02ad21d70285e3cec7204f6d8933fa07fdde16").unwrap(),
            from_hex("0x025f84c58f3d7cdc2f8b021f130aee1b89d6edcb43139f2cca0fa6ae8caa86ff").unwrap(),
            from_hex("0x00cf44cb881d539fba196c90a629c159b57974850d8e650b765efbdfb014101c").unwrap(),
        ],
        vec![
            from_hex("0x0396dbdc13a64c52668aacf6409f8bc6abba447e15dcb5690e7b670602bb9e91").unwrap(),
            from_hex("0x010ffe7154727122d7026d92e8d4e585adc9807b62310883c0c6502ffe2cd087").unwrap(),
            from_hex("0x03149e4fb149e713a5360f0fb97e8a512d0fc86ae806a40acae092bd7c69900d").unwrap(),
        ],
        vec![
            from_hex("0x02bb6601bdbec2fe689fddabb59a645ee08de0e545a6d215fdc52620ad72eeb5").unwrap(),
            from_hex("0x011e5eff841478933520b03068fd84769ce1446ab9b5bfa6b2a75c33fb4f58a8").unwrap(),
            from_hex("0x024e5dfbda4e272e991e4b02f30889baa5d0c5d49bab2c5ad1c063c60c80a6c5").unwrap(),
        ],
        vec![
            from_hex("0x00aa5c75ee6bd14a7f41be8a8c385b45c0ae9e1055188f73b17d347c12d39260").unwrap(),
            from_hex("0x001cc4d0d20518d1221dd44ba295b397689c6eecf189eef70120dcb10da30bb4").unwrap(),
            from_hex("0x026bf5c30b55e038ba1a238b3bf9355f4d401371d20cd798cb0c116a2c05a10b").unwrap(),
        ],
        vec![
            from_hex("0x03677c11de06dff458e2e75a9d2f3fb0d6fd3c48eec74034ce7c1f2b2b0927a4").unwrap(),
            from_hex("0x02aa6924ad56a5824bb3875e6a8c47476dc81e3952e3eff025ee1b84f7f0c381").unwrap(),
            from_hex("0x01a3032aff572f3baf6f0d7e62f986af5d2e64231af84e25967f8ae5c225c213").unwrap(),
        ],
        vec![
            from_hex("0x01c80ef5070b91391a2eb8ee0ce94f996468a33f5ada53cd30c86bad718fbc25").unwrap(),
            from_hex("0x0353fed6ac39e84209acb97e8819d7045643d2c724fba3ffe292803eda97a0c0").unwrap(),
            from_hex("0x03017359916d40f9e6575af6d35d2b190649cceb26c48907a972bfec67796ee1").unwrap(),
        ],
        vec![
            from_hex("0x0055ebd0023afb6f8ba8d8fdca879148fa08e1e776884b5f11cf458e74fc6d29").unwrap(),
            from_hex("0x00d791bb67a3273b06150afd059cce00779483a27f84521527a10b55b0c59996").unwrap(),
            from_hex("0x00e69c9d83dfb0ae81ddfe054f0d5fb4f33b762282dcb24bc97261d3946d4ad7").unwrap(),
        ],
        vec![
            from_hex("0x0010efc2de030d91686881d9d9080ae757d76d573cf629e5ba747b3b39f5486e").unwrap(),
            from_hex("0x0118deda5ada35552248404952bc45005423eef890074272a49d6b848310cae7").unwrap(),
            from_hex("0x03d3e5dd62361d59ac54fa05121a14b36fa6cc7ae9e363513c90a141c1e57061").unwrap(),
        ],
        vec![
            from_hex("0x01dda8937b71763234e1aa685a5b7361c1ab1e3d30d50f18f10986d7241f586d").unwrap(),
            from_hex("0x03690d5d0cd0689fb65aaa117fab941b13061a489c7ad3ce9229175384259689").unwrap(),
            from_hex("0x014fbf9d4f3f10454e9248dae27fc005482ad1e01986ff446ff2c9d1a2bb6a85").unwrap(),
        ],
        vec![
            from_hex("0x0090a06107fa03f9d00476d672a7456ffd388b64699051db7dee37ff4c7e6767").unwrap(),
            from_hex("0x01a1deacae3634393bae2ed06ebf283376ae68ea8089ca0c9edf49b48c2e6bae").unwrap(),
            from_hex("0x00df2871cc32d0e73049c1343feb396076f1d357916763972efd07e3292ed512").unwrap(),
        ],
        vec![
            from_hex("0x008e6c962fd2bec44aaafabf4aac16eb7715fde75f45306b539b51b6823bc315").unwrap(),
            from_hex("0x02bb848c4a36c1cdef6b15b181942cb24083ad68278ec7ccf480756d175d7d9a").unwrap(),
            from_hex("0x01a2077d541ca1431ee51331cde405337558021dfa7425d9169b218ebd7dd1bb").unwrap(),
        ],
        vec![
            from_hex("0x029b65931beda4cedff9d97405915bcbc6bf6c92052b6444b06bcd0f1682a8cf").unwrap(),
            from_hex("0x01a495989061804d2bf82a1b85a6bad80cd8c196116086566c0bb7e7dfbc944b").unwrap(),
            from_hex("0x03bd7178b20e3a8ea4b871981cfab8dc85761028bfe7641a1d528d229f6fae38").unwrap(),
        ],
        vec![
            from_hex("0x03e644569e23fba86891f59caffdff53fc3f0b420a6909b63e66a78a57f96c65").unwrap(),
            from_hex("0x033a94e83fced654deb817c02a898d826cf65c68712a6577dba1611f7bde436b").unwrap(),
            from_hex("0x0040630431ec89f00b53d83bb9f5cb8a41e9fc9bf5e46acf749b5bfd3d2fae70").unwrap(),
        ],
        vec![
            from_hex("0x03b9f636a08b9c5d3f1fe4b58b2eb5ebdcf33d04c330229ecaf42553f42f4d40").unwrap(),
            from_hex("0x02d7fbfe825cbf6fd6ea6ba665b167608a22913971715b055c60b71c54f350dc").unwrap(),
            from_hex("0x00771b6c67a2454a32c73fcc76e5259fb7c933eb82a5214dd3969b997a6a4652").unwrap(),
        ],
        vec![
            from_hex("0x00431d92909da2f0b71f03964c597733b24a14a4b817d07ee4ebc641f88bf9fc").unwrap(),
            from_hex("0x03ebafd4f90a014580406535e259bead22018c1f099efdbb357e1d9b17b6c2a2").unwrap(),
            from_hex("0x016664cbc87543f93d50918cf06511cd6d79039f4c2346ffa6a5b6a1c521d846").unwrap(),
        ],
        vec![
            from_hex("0x014c4bb5ed7c1819d5ae71b929b499e9b689989e8b0192a352b9365bc5703f94").unwrap(),
            from_hex("0x035b9d7f7a1da1d87bbb370139ec40f32c6e142a544199c9e9fcf0e0d2da5594").unwrap(),
            from_hex("0x008d5fb7e7ec4655e1a06d023c1af89a471762219ec2989aef58d0a7c160514d").unwrap(),
        ],
        vec![
            from_hex("0x0052543df0b5b33ebab3c20ac8d97e1b82ed64f403c01bcadd81b9c9a3753f0f").unwrap(),
            from_hex("0x00047520c96ea10a4dce4b8badca9528bbeeba427fafe69a892056bbbc7c7a72").unwrap(),
            from_hex("0x00eebb53a17d362b02492d798684051a76859a20fb3e58146d4abb63cc1041c2").unwrap(),
        ],
        vec![
            from_hex("0x00791cbe16bbc4ae7fa4b7b02d9fc8ef98d635cd3e9728afe3aeab23139377a9").unwrap(),
            from_hex("0x01e32c4f6012453d23f64e9210f6f9c5c97b671f87dddd58d414ed68b3c9a21a").unwrap(),
            from_hex("0x0129368c89b5e46e83bbc6be1e4adee56f0ccfcc000bd101d1c0b518cf9cbc4f").unwrap(),
        ],
        vec![
            from_hex("0x03eb7b0f647797488c1dd1e50180816d356d138e18d68ec7ff91c41ce6f4884e").unwrap(),
            from_hex("0x013cb99029d915f3d74deb42f7cf404a6ebdcd0f2f46dda2ebf76b43f8dfd749").unwrap(),
            from_hex("0x03c5033b0f286eb72a1193535cee97b8469a3ada8845f6407400fb2faa614f17").unwrap(),
        ],
        vec![
            from_hex("0x0050161b27849ee9bfb5cd729f2f99a087190ba78e2feadfec04e01cf7b3daef").unwrap(),
            from_hex("0x01daf4e4171e4ea3d115e9b8ac8b8e564ae3daf1e50fd0dd33da4db9f9598b57").unwrap(),
            from_hex("0x019802327513cb48230a597ccd789327bbe47eda84e630c71446874a382d2ccc").unwrap(),
        ],
        vec![
            from_hex("0x0144efa46e25fca6c30b8408d46dc0300e74e9a5e18554b188e3bc35aaab1e1d").unwrap(),
            from_hex("0x02945618652759f444679d9dbee7e7b4c7dbffffe642abe1781044a9a2f4424e").unwrap(),
            from_hex("0x01e4dae0c91578c7ba9cc9e45aa350dfe8690ca6e5e906352d9c9b93c81c47f0").unwrap(),
        ],
        vec![
            from_hex("0x01f6968b85035a2fff54b67d3e39d09d69cc541ced158d9f0a413601c247d095").unwrap(),
            from_hex("0x00b5a43934f200459b22dcc59d8ac131a81670568d3f601793ac77a9b8ed9ca8").unwrap(),
            from_hex("0x03daa74a9faadb70f10d720f77ffe50a00657c2abab80d9793b4f0b745adaba6").unwrap(),
        ],
        vec![
            from_hex("0x02bc6f5e10124e45c649cd0ea2cc789e81a969a32ee62ceffb4f6e923c747012").unwrap(),
            from_hex("0x01edb3e77c647bf470490c9e3fce3e8a6f029396c4d6f5c9f85633de235251f6").unwrap(),
            from_hex("0x03bd4dec479b0055916b557c3e9649de9db61a97773c574927a8a7b09fac57df").unwrap(),
        ],
        vec![
            from_hex("0x031421777c106c76d0701c032dea9e3081d953d5c46af9f5df5dfb0336520788").unwrap(),
            from_hex("0x008390c462ee68d2cc16dbe096e7d820563171ec447f6df6164b922dc43399c2").unwrap(),
            from_hex("0x00938401840129b7fb50a9fcb71611ad5607e40064908238e05b78c91f74b67e").unwrap(),
        ],
        vec![
            from_hex("0x00baa2ad97c86efa5e5f537aac34a8ed920f5fd69fcac58eed9ed6e15d75ae0c").unwrap(),
            from_hex("0x03559803c902ac640a8d0f606ac4e6d7a897bbe9f4005f6c80941207e7b22484").unwrap(),
            from_hex("0x03df547e920ccb9a981534b030cd2ccf7135b3602f960daad2e4c91709f1ad2f").unwrap(),
        ],
        vec![
            from_hex("0x026bd857d965b19ed7669eef9bc69e8309fa29f32e854db7a2d8dec74595f6dd").unwrap(),
            from_hex("0x02af7f0ddf81fec89c7b89a18c8f2bceae0d52f345919d3bb7e2dcdc545bad86").unwrap(),
            from_hex("0x0069678f69c4031ecda42fc7ed0ede5cbf4ca9f9dc22e423cd8c26d26d46e62c").unwrap(),
        ],
        vec![
            from_hex("0x0279a21d80075cc540f313a103f1c793b4f8b2d4329f1e4f52360666739769c1").unwrap(),
            from_hex("0x01ea8c0a796419033424db8909b22baf5b894e03b6697159ae128bb4441071d2").unwrap(),
            from_hex("0x0366b003f5affc4664e1f32e8a332b4257334ef08c9777249396a60b437abd53").unwrap(),
        ],
        vec![
            from_hex("0x03f38ea496ff522369938f662b64cb25b818e8290c9469f47f519b237ec63180").unwrap(),
            from_hex("0x00e2fff5e11013974485e03f052168ed24d7f3cf5ba938591e445e0526a77370").unwrap(),
            from_hex("0x037546f93699b435ceb61c624666068858fc356e9818371bda061924fd2954c6").unwrap(),
        ],
        vec![
            from_hex("0x017d19c6790e1aaf25ca5c9f8286527b628d6cd3874c6b633df4be080570710b").unwrap(),
            from_hex("0x004004525f018a4812911762b2e2d84fe9a722e417ac7a8f4e53042eae306716").unwrap(),
            from_hex("0x03cf92b5c8c780a76c13e967b181895e160f3ce3e97cc17ef7500ee04d3f7460").unwrap(),
        ],
    ];
    pub static ref RESCUE_ST_PARAMS: Arc<RescueParams<Scalar>> = Arc::new(RescueParams::new(
        3,
        3,
        [
            0xaaaaaaaaaaaaaaab,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0x2a6aaaaaaaaaaaa,
        ],
        22,
        &MDS3,
        &RC3
    ));
}
