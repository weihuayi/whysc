
namespace WHYSC{
namespace Quadrature{

const double TetrahedronQuadratureData[210][5] ={

        /**
         * \brife 一阶积分
         */
        {0.2500000000000000,	0.2500000000000000,	0.2500000000000000,	0.2500000000000000,	1.0000000000000000},

        /**
         * \brife 二阶积分
         */
        {0.5854101966249680,	0.1381966011250110,	0.1381966011250110,	0.1381966011250110,	0.2500000000000000},
        {0.1381966011250110,	0.5854101966249680,	0.1381966011250110,	0.1381966011250110,	0.2500000000000000},
        {0.1381966011250110,	0.1381966011250110,	0.5854101966249680,	0.1381966011250110,	0.2500000000000000},
        {0.1381966011250110,	0.1381966011250110,	0.1381966011250110,	0.5854101966249680,	0.2500000000000000},

        /**
         * \brife 三阶积分
         */
        {0.7784952948213300,	0.0738349017262234,	0.0738349017262234,	0.0738349017262234,	0.0476331348432089},
        {0.0738349017262234,	0.7784952948213300,	0.0738349017262234,	0.0738349017262234,	0.0476331348432089},
        {0.0738349017262234,	0.0738349017262234,	0.7784952948213300,	0.0738349017262234,	0.0476331348432089},
        {0.0738349017262234,	0.0738349017262234,	0.0738349017262234,	0.7784952948213300,	0.0476331348432089},
        {0.4062443438840510,	0.4062443438840510,	0.0937556561159491,	0.0937556561159491,	0.1349112434378610},
        {0.4062443438840510,	0.0937556561159491,	0.4062443438840510,	0.0937556561159491,	0.1349112434378610},
        {0.4062443438840510,	0.0937556561159491,	0.0937556561159491,	0.4062443438840510,	0.1349112434378610},
        {0.0937556561159491,	0.4062443438840510,	0.4062443438840510,	0.0937556561159491,	0.1349112434378610},
        {0.0937556561159491,	0.4062443438840510,	0.0937556561159491,	0.4062443438840510,	0.1349112434378610},
        {0.0937556561159491,	0.0937556561159491,	0.4062443438840510,	0.4062443438840510,	0.1349112434378610},

        /**
         * \brife 四阶积分
         */
        {0.9029422158182680,	0.0323525947272439,	0.0323525947272439,	0.0323525947272439,	0.0070670747944695},
        {0.0323525947272439,	0.9029422158182680,	0.0323525947272439,	0.0323525947272439,	0.0070670747944695},
        {0.0323525947272439,	0.0323525947272439,	0.9029422158182680,	0.0323525947272439,	0.0070670747944695},
        {0.0323525947272439,	0.0323525947272439,	0.0323525947272439,	0.9029422158182680,	0.0070670747944695},
        {0.2626825838877790,	0.6165965330619370,	0.0603604415251421,	0.0603604415251421,	0.0469986689718877},
        {0.6165965330619370,	0.2626825838877790,	0.0603604415251421,	0.0603604415251421,	0.0469986689718877},
        {0.2626825838877790,	0.0603604415251421,	0.6165965330619370,	0.0603604415251421,	0.0469986689718877},
        {0.6165965330619370,	0.0603604415251421,	0.2626825838877790,	0.0603604415251421,	0.0469986689718877},
        {0.2626825838877790,	0.0603604415251421,	0.0603604415251421,	0.6165965330619370,	0.0469986689718877},
        {0.6165965330619370,	0.0603604415251421,	0.0603604415251421,	0.2626825838877790,	0.0469986689718877},
        {0.0603604415251421,	0.2626825838877790,	0.6165965330619370,	0.0603604415251421,	0.0469986689718877},
        {0.0603604415251421,	0.6165965330619370,	0.2626825838877790,	0.0603604415251421,	0.0469986689718877},
        {0.0603604415251421,	0.2626825838877790,	0.0603604415251421,	0.6165965330619370,	0.0469986689718877},
        {0.0603604415251421,	0.6165965330619370,	0.0603604415251421,	0.2626825838877790,	0.0469986689718877},
        {0.0603604415251421,	0.0603604415251421,	0.2626825838877790,	0.6165965330619370,	0.0469986689718877},
        {0.0603604415251421,	0.0603604415251421,	0.6165965330619370,	0.2626825838877790,	0.0469986689718877},
        {0.3097693042728620,	0.3097693042728620,	0.3097693042728620,	0.0706920871814129,	0.1019369182898680},
        {0.3097693042728620,	0.3097693042728620,	0.0706920871814129,	0.3097693042728620,	0.1019369182898680},
        {0.3097693042728620,	0.0706920871814129,	0.3097693042728620,	0.3097693042728620,	0.1019369182898680},
        {0.0706920871814129,	0.3097693042728620,	0.3097693042728620,	0.3097693042728620,	0.1019369182898680},

        /**
         * \brife 五阶积分
         */
        {0.9197896733368800,	0.0267367755543735,	0.0267367755543735,	0.0267367755543735,	0.0021900463965388},
        {0.0267367755543735,	0.9197896733368800,	0.0267367755543735,	0.0267367755543735,	0.0021900463965388},
        {0.0267367755543735,	0.0267367755543735,	0.9197896733368800,	0.0267367755543735,	0.0021900463965388},
        {0.0267367755543735,	0.0267367755543735,	0.0267367755543735,	0.9197896733368800,	0.0021900463965388},
        {0.1740356302468940,	0.7477598884818090,	0.0391022406356488,	0.0391022406356488,	0.0143395670177665},
        {0.7477598884818090,	0.1740356302468940,	0.0391022406356488,	0.0391022406356488,	0.0143395670177665},
        {0.1740356302468940,	0.0391022406356488,	0.7477598884818090,	0.0391022406356488,	0.0143395670177665},
        {0.7477598884818090,	0.0391022406356488,	0.1740356302468940,	0.0391022406356488,	0.0143395670177665},
        {0.1740356302468940,	0.0391022406356488,	0.0391022406356488,	0.7477598884818090,	0.0143395670177665},
        {0.7477598884818090,	0.0391022406356488,	0.0391022406356488,	0.1740356302468940,	0.0143395670177665},
        {0.0391022406356488,	0.1740356302468940,	0.7477598884818090,	0.0391022406356488,	0.0143395670177665},
        {0.0391022406356488,	0.7477598884818090,	0.1740356302468940,	0.0391022406356488,	0.0143395670177665},
        {0.0391022406356488,	0.1740356302468940,	0.0391022406356488,	0.7477598884818090,	0.0143395670177665},
        {0.0391022406356488,	0.7477598884818090,	0.0391022406356488,	0.1740356302468940,	0.0143395670177665},
        {0.0391022406356488,	0.0391022406356488,	0.1740356302468940,	0.7477598884818090,	0.0143395670177665},
        {0.0391022406356488,	0.0391022406356488,	0.7477598884818090,	0.1740356302468940,	0.0143395670177665},
        {0.4547545999844830,	0.4547545999844830,	0.0452454000155172,	0.0452454000155172,	0.0250305395686746},
        {0.4547545999844830,	0.0452454000155172,	0.4547545999844830,	0.0452454000155172,	0.0250305395686746},
        {0.4547545999844830,	0.0452454000155172,	0.0452454000155172,	0.4547545999844830,	0.0250305395686746},
        {0.0452454000155172,	0.4547545999844830,	0.4547545999844830,	0.0452454000155172,	0.0250305395686746},
        {0.0452454000155172,	0.4547545999844830,	0.0452454000155172,	0.4547545999844830,	0.0250305395686746},
        {0.0452454000155172,	0.0452454000155172,	0.4547545999844830,	0.4547545999844830,	0.0250305395686746},
        {0.5031186450145980,	0.2232010379623150,	0.2232010379623150,	0.0504792790607720,	0.0479839333057554},
        {0.2232010379623150,	0.5031186450145980,	0.2232010379623150,	0.0504792790607720,	0.0479839333057554},
        {0.2232010379623150,	0.2232010379623150,	0.5031186450145980,	0.0504792790607720,	0.0479839333057554},
        {0.5031186450145980,	0.2232010379623150,	0.0504792790607720,	0.2232010379623150,	0.0479839333057554},
        {0.2232010379623150,	0.5031186450145980,	0.0504792790607720,	0.2232010379623150,	0.0479839333057554},
        {0.2232010379623150,	0.2232010379623150,	0.0504792790607720,	0.5031186450145980,	0.0479839333057554},
        {0.5031186450145980,	0.0504792790607720,	0.2232010379623150,	0.2232010379623150,	0.0479839333057554},
        {0.2232010379623150,	0.0504792790607720,	0.5031186450145980,	0.2232010379623150,	0.0479839333057554},
        {0.2232010379623150,	0.0504792790607720,	0.2232010379623150,	0.5031186450145980,	0.0479839333057554},
        {0.0504792790607720,	0.5031186450145980,	0.2232010379623150,	0.2232010379623150,	0.0479839333057554},
        {0.0504792790607720,	0.2232010379623150,	0.5031186450145980,	0.2232010379623150,	0.0479839333057554},
        {0.0504792790607720,	0.2232010379623150,	0.2232010379623150,	0.5031186450145980,	0.0479839333057554},
        {0.2500000000000000,	0.2500000000000000,	0.2500000000000000,	0.2500000000000000,	0.0931745731195340},

        /**
         * \brife 六阶积分
         */
        {0.9551438045408220,	0.0149520651530592,	0.0149520651530592,	0.0149520651530592,	0.0010373112336140},
        {0.0149520651530592,	0.9551438045408220,	0.0149520651530592,	0.0149520651530592,	0.0010373112336140},
        {0.0149520651530592,	0.0149520651530592,	0.9551438045408220,	0.0149520651530592,	0.0010373112336140},
        {0.0149520651530592,	0.0149520651530592,	0.0149520651530592,	0.9551438045408220,	0.0010373112336140},
        {0.7799760084415400,	0.1518319491659370,	0.0340960211962615,	0.0340960211962615,	0.0096016645399480},
        {0.1518319491659370,	0.7799760084415400,	0.0340960211962615,	0.0340960211962615,	0.0096016645399480},
        {0.7799760084415400,	0.0340960211962615,	0.1518319491659370,	0.0340960211962615,	0.0096016645399480},
        {0.1518319491659370,	0.0340960211962615,	0.7799760084415400,	0.0340960211962615,	0.0096016645399480},
        {0.7799760084415400,	0.0340960211962615,	0.0340960211962615,	0.1518319491659370,	0.0096016645399480},
        {0.1518319491659370,	0.0340960211962615,	0.0340960211962615,	0.7799760084415400,	0.0096016645399480},
        {0.0340960211962615,	0.7799760084415400,	0.1518319491659370,	0.0340960211962615,	0.0096016645399480},
        {0.0340960211962615,	0.1518319491659370,	0.7799760084415400,	0.0340960211962615,	0.0096016645399480},
        {0.0340960211962615,	0.7799760084415400,	0.0340960211962615,	0.1518319491659370,	0.0096016645399480},
        {0.0340960211962615,	0.1518319491659370,	0.0340960211962615,	0.7799760084415400,	0.0096016645399480},
        {0.0340960211962615,	0.0340960211962615,	0.7799760084415400,	0.1518319491659370,	0.0096016645399480},
        {0.0340960211962615,	0.0340960211962615,	0.1518319491659370,	0.7799760084415400,	0.0096016645399480},
        {0.3549340560639790,	0.5526556431060170,	0.0462051504150017,	0.0462051504150017,	0.0164493976798232},
        {0.5526556431060170,	0.3549340560639790,	0.0462051504150017,	0.0462051504150017,	0.0164493976798232},
        {0.3549340560639790,	0.0462051504150017,	0.5526556431060170,	0.0462051504150017,	0.0164493976798232},
        {0.5526556431060170,	0.0462051504150017,	0.3549340560639790,	0.0462051504150017,	0.0164493976798232},
        {0.3549340560639790,	0.0462051504150017,	0.0462051504150017,	0.5526556431060170,	0.0164493976798232},
        {0.5526556431060170,	0.0462051504150017,	0.0462051504150017,	0.3549340560639790,	0.0164493976798232},
        {0.0462051504150017,	0.3549340560639790,	0.5526556431060170,	0.0462051504150017,	0.0164493976798232},
        {0.0462051504150017,	0.5526556431060170,	0.3549340560639790,	0.0462051504150017,	0.0164493976798232},
        {0.0462051504150017,	0.3549340560639790,	0.0462051504150017,	0.5526556431060170,	0.0164493976798232},
        {0.0462051504150017,	0.5526556431060170,	0.0462051504150017,	0.3549340560639790,	0.0164493976798232},
        {0.0462051504150017,	0.0462051504150017,	0.3549340560639790,	0.5526556431060170,	0.0164493976798232},
        {0.0462051504150017,	0.0462051504150017,	0.5526556431060170,	0.3549340560639790,	0.0164493976798232},
        {0.5381043228880020,	0.2281904610687610,	0.2281904610687610,	0.0055147549744775,	0.0153747766513310},
        {0.2281904610687610,	0.5381043228880020,	0.2281904610687610,	0.0055147549744775,	0.0153747766513310},
        {0.2281904610687610,	0.2281904610687610,	0.5381043228880020,	0.0055147549744775,	0.0153747766513310},
        {0.5381043228880020,	0.2281904610687610,	0.0055147549744775,	0.2281904610687610,	0.0153747766513310},
        {0.2281904610687610,	0.5381043228880020,	0.0055147549744775,	0.2281904610687610,	0.0153747766513310},
        {0.2281904610687610,	0.2281904610687610,	0.0055147549744775,	0.5381043228880020,	0.0153747766513310},
        {0.5381043228880020,	0.0055147549744775,	0.2281904610687610,	0.2281904610687610,	0.0153747766513310},
        {0.2281904610687610,	0.0055147549744775,	0.5381043228880020,	0.2281904610687610,	0.0153747766513310},
        {0.2281904610687610,	0.0055147549744775,	0.2281904610687610,	0.5381043228880020,	0.0153747766513310},
        {0.0055147549744775,	0.5381043228880020,	0.2281904610687610,	0.2281904610687610,	0.0153747766513310},
        {0.0055147549744775,	0.2281904610687610,	0.5381043228880020,	0.2281904610687610,	0.0153747766513310},
        {0.0055147549744775,	0.2281904610687610,	0.2281904610687610,	0.5381043228880020,	0.0153747766513310},
        {0.1961837595745600,	0.3523052600879940,	0.3523052600879940,	0.0992057202494530,	0.0293520118375230},
        {0.3523052600879940,	0.1961837595745600,	0.3523052600879940,	0.0992057202494530,	0.0293520118375230},
        {0.3523052600879940,	0.3523052600879940,	0.1961837595745600,	0.0992057202494530,	0.0293520118375230},
        {0.1961837595745600,	0.3523052600879940,	0.0992057202494530,	0.3523052600879940,	0.0293520118375230},
        {0.3523052600879940,	0.1961837595745600,	0.0992057202494530,	0.3523052600879940,	0.0293520118375230},
        {0.3523052600879940,	0.3523052600879940,	0.0992057202494530,	0.1961837595745600,	0.0293520118375230},
        {0.1961837595745600,	0.0992057202494530,	0.3523052600879940,	0.3523052600879940,	0.0293520118375230},
        {0.3523052600879940,	0.0992057202494530,	0.1961837595745600,	0.3523052600879940,	0.0293520118375230},
        {0.3523052600879940,	0.0992057202494530,	0.3523052600879940,	0.1961837595745600,	0.0293520118375230},
        {0.0992057202494530,	0.1961837595745600,	0.3523052600879940,	0.3523052600879940,	0.0293520118375230},
        {0.0992057202494530,	0.3523052600879940,	0.1961837595745600,	0.3523052600879940,	0.0293520118375230},
        {0.0992057202494530,	0.3523052600879940,	0.3523052600879940,	0.1961837595745600,	0.0293520118375230},
        {0.5965649956210170,	0.1344783347929940,	0.1344783347929940,	0.1344783347929940,	0.0366291366405108},
        {0.1344783347929940,	0.5965649956210170,	0.1344783347929940,	0.1344783347929940,	0.0366291366405108},
        {0.1344783347929940,	0.1344783347929940,	0.5965649956210170,	0.1344783347929940,	0.0366291366405108},
        {0.1344783347929940,	0.1344783347929940,	0.1344783347929940,	0.5965649956210170,	0.0366291366405108},

        /**
         * \brife 七阶积分
         */
        {0.9193645767555490,	0.0268784744148170,	0.0268784744148170,	0.0268784744148170,	0.0021449351443160},
        {0.0268784744148170,	0.9193645767555490,	0.0268784744148170,	0.0268784744148170,	0.0021449351443160},
        {0.0268784744148170,	0.0268784744148170,	0.9193645767555490,	0.0268784744148170,	0.0021449351443160},
        {0.0268784744148170,	0.0268784744148170,	0.0268784744148170,	0.9193645767555490,	0.0021449351443160},
        {0.4385779725895910,	0.1871406758034700,	0.1871406758034700,	0.1871406758034700,	0.0208266416907690},
        {0.1871406758034700,	0.4385779725895910,	0.1871406758034700,	0.1871406758034700,	0.0208266416907690},
        {0.1871406758034700,	0.1871406758034700,	0.4385779725895910,	0.1871406758034700,	0.0208266416907690},
        {0.1871406758034700,	0.1871406758034700,	0.1871406758034700,	0.4385779725895910,	0.0208266416907690},
        {0.4735758351279370,	0.4735758351279370,	0.0264241648720630,	0.0264241648720630,	0.0072101360644550},
        {0.4735758351279370,	0.0264241648720630,	0.4735758351279370,	0.0264241648720630,	0.0072101360644550},
        {0.4735758351279370,	0.0264241648720630,	0.0264241648720630,	0.4735758351279370,	0.0072101360644550},
        {0.0264241648720630,	0.4735758351279370,	0.4735758351279370,	0.0264241648720630,	0.0072101360644550},
        {0.0264241648720630,	0.4735758351279370,	0.0264241648720630,	0.4735758351279370,	0.0072101360644550},
        {0.0264241648720630,	0.0264241648720630,	0.4735758351279370,	0.4735758351279370,	0.0072101360644550},
        {0.3520452620273560,	0.3520452620273560,	0.1479547379726440,	0.1479547379726440,	0.0307989191597120},
        {0.3520452620273560,	0.1479547379726440,	0.3520452620273560,	0.1479547379726440,	0.0307989191597120},
        {0.3520452620273560,	0.1479547379726440,	0.1479547379726440,	0.3520452620273560,	0.0307989191597120},
        {0.1479547379726440,	0.3520452620273560,	0.3520452620273560,	0.1479547379726440,	0.0307989191597120},
        {0.1479547379726440,	0.3520452620273560,	0.1479547379726440,	0.3520452620273560,	0.0307989191597120},
        {0.1479547379726440,	0.1479547379726440,	0.3520452620273560,	0.3520452620273560,	0.0307989191597120},
        {0.2257832058669400,	0.7323099096929470,	0.0209534422200560,	0.0209534422200560,	0.0043578448138640},
        {0.2257832058669400,	0.0209534422200560,	0.7323099096929470,	0.0209534422200560,	0.0043578448138640},
        {0.2257832058669400,	0.0209534422200560,	0.0209534422200560,	0.7323099096929470,	0.0043578448138640},
        {0.0209534422200560,	0.2257832058669400,	0.7323099096929470,	0.0209534422200560,	0.0043578448138640},
        {0.0209534422200560,	0.2257832058669400,	0.0209534422200560,	0.7323099096929470,	0.0043578448138640},
        {0.0209534422200560,	0.0209534422200560,	0.2257832058669400,	0.7323099096929470,	0.0043578448138640},
        {0.7323099096929470,	0.2257832058669400,	0.0209534422200560,	0.0209534422200560,	0.0043578448138640},
        {0.7323099096929470,	0.0209534422200560,	0.2257832058669400,	0.0209534422200560,	0.0043578448138640},
        {0.7323099096929470,	0.0209534422200560,	0.0209534422200560,	0.2257832058669400,	0.0043578448138640},
        {0.0209534422200560,	0.7323099096929470,	0.2257832058669400,	0.0209534422200560,	0.0043578448138640},
        {0.0209534422200560,	0.7323099096929470,	0.0209534422200560,	0.2257832058669400,	0.0043578448138640},
        {0.0209534422200560,	0.0209534422200560,	0.7323099096929470,	0.2257832058669400,	0.0043578448138640},
        {0.1584629396660920,	0.6475575940869750,	0.0969897331234660,	0.0969897331234660,	0.0085935306778330},
        {0.1584629396660920,	0.0969897331234660,	0.6475575940869750,	0.0969897331234660,	0.0085935306778330},
        {0.1584629396660920,	0.0969897331234660,	0.0969897331234660,	0.6475575940869750,	0.0085935306778330},
        {0.0969897331234660,	0.1584629396660920,	0.6475575940869750,	0.0969897331234660,	0.0085935306778330},
        {0.0969897331234660,	0.1584629396660920,	0.0969897331234660,	0.6475575940869750,	0.0085935306778330},
        {0.0969897331234660,	0.0969897331234660,	0.1584629396660920,	0.6475575940869750,	0.0085935306778330},
        {0.6475575940869750,	0.1584629396660920,	0.0969897331234660,	0.0969897331234660,	0.0085935306778330},
        {0.6475575940869750,	0.0969897331234660,	0.1584629396660920,	0.0969897331234660,	0.0085935306778330},
        {0.6475575940869750,	0.0969897331234660,	0.0969897331234660,	0.1584629396660920,	0.0085935306778330},
        {0.0969897331234660,	0.6475575940869750,	0.1584629396660920,	0.0969897331234660,	0.0085935306778330},
        {0.0969897331234660,	0.6475575940869750,	0.0969897331234660,	0.1584629396660920,	0.0085935306778330},
        {0.0969897331234660,	0.0969897331234660,	0.6475575940869750,	0.1584629396660920,	0.0085935306778330},
        {0.3221114318308570,	0.3221114318308570,	0.3221114318308570,	0.0336657045074290,	0.0230006816692860},
        {0.3221114318308570,	0.3221114318308570,	0.0336657045074290,	0.3221114318308570,	0.0230006816692860},
        {0.3221114318308570,	0.0336657045074290,	0.3221114318308570,	0.3221114318308570,	0.0230006816692860},
        {0.0336657045074290,	0.3221114318308570,	0.3221114318308570,	0.3221114318308570,	0.0230006816692860},
        {0.7929392564696180,	0.0976081628904420,	0.0976081628904420,	0.0118444177494980,	0.0048630639049120},
        {0.0976081628904420,	0.7929392564696180,	0.0976081628904420,	0.0118444177494980,	0.0048630639049120},
        {0.0976081628904420,	0.0976081628904420,	0.7929392564696180,	0.0118444177494980,	0.0048630639049120},
        {0.7929392564696180,	0.0976081628904420,	0.0118444177494980,	0.0976081628904420,	0.0048630639049120},
        {0.0976081628904420,	0.7929392564696180,	0.0118444177494980,	0.0976081628904420,	0.0048630639049120},
        {0.0976081628904420,	0.0976081628904420,	0.0118444177494980,	0.7929392564696180,	0.0048630639049120},
        {0.7929392564696180,	0.0118444177494980,	0.0976081628904420,	0.0976081628904420,	0.0048630639049120},
        {0.0976081628904420,	0.0118444177494980,	0.7929392564696180,	0.0976081628904420,	0.0048630639049120},
        {0.0976081628904420,	0.0118444177494980,	0.0976081628904420,	0.7929392564696180,	0.0048630639049120},
        {0.0118444177494980,	0.7929392564696180,	0.0976081628904420,	0.0976081628904420,	0.0048630639049120},
        {0.0118444177494980,	0.0976081628904420,	0.7929392564696180,	0.0976081628904420,	0.0048630639049120},
        {0.0118444177494980,	0.0976081628904420,	0.0976081628904420,	0.7929392564696180,	0.0048630639049120},
        {0.5411844128002370,	0.1335581607035680,	0.2965010205431240,	0.0287564059530710,	0.0155951400782590},
        {0.2965010205431240,	0.5411844128002370,	0.1335581607035680,	0.0287564059530710,	0.0155951400782590},
        {0.1335581607035680,	0.2965010205431240,	0.5411844128002370,	0.0287564059530710,	0.0155951400782590},
        {0.5411844128002370,	0.2965010205431240,	0.1335581607035680,	0.0287564059530710,	0.0155951400782590},
        {0.1335581607035680,	0.5411844128002370,	0.2965010205431240,	0.0287564059530710,	0.0155951400782590},
        {0.2965010205431240,	0.1335581607035680,	0.5411844128002370,	0.0287564059530710,	0.0155951400782590},
        {0.5411844128002370,	0.1335581607035680,	0.0287564059530710,	0.2965010205431240,	0.0155951400782590},
        {0.2965010205431240,	0.5411844128002370,	0.0287564059530710,	0.1335581607035680,	0.0155951400782590},
        {0.1335581607035680,	0.2965010205431240,	0.0287564059530710,	0.5411844128002370,	0.0155951400782590},
        {0.5411844128002370,	0.2965010205431240,	0.0287564059530710,	0.1335581607035680,	0.0155951400782590},
        {0.1335581607035680,	0.5411844128002370,	0.0287564059530710,	0.2965010205431240,	0.0155951400782590},
        {0.2965010205431240,	0.1335581607035680,	0.0287564059530710,	0.5411844128002370,	0.0155951400782590},
        {0.5411844128002370,	0.0287564059530710,	0.1335581607035680,	0.2965010205431240,	0.0155951400782590},
        {0.2965010205431240,	0.0287564059530710,	0.5411844128002370,	0.1335581607035680,	0.0155951400782590},
        {0.1335581607035680,	0.0287564059530710,	0.2965010205431240,	0.5411844128002370,	0.0155951400782590},
        {0.5411844128002370,	0.0287564059530710,	0.2965010205431240,	0.1335581607035680,	0.0155951400782590},
        {0.1335581607035680,	0.0287564059530710,	0.5411844128002370,	0.2965010205431240,	0.0155951400782590},
        {0.2965010205431240,	0.0287564059530710,	0.1335581607035680,	0.5411844128002370,	0.0155951400782590},
        {0.0287564059530710,	0.5411844128002370,	0.1335581607035680,	0.2965010205431240,	0.0155951400782590},
        {0.0287564059530710,	0.2965010205431240,	0.5411844128002370,	0.1335581607035680,	0.0155951400782590},
        {0.0287564059530710,	0.1335581607035680,	0.2965010205431240,	0.5411844128002370,	0.0155951400782590},
        {0.0287564059530710,	0.5411844128002370,	0.2965010205431240,	0.1335581607035680,	0.0155951400782590},
        {0.0287564059530710,	0.1335581607035680,	0.5411844128002370,	0.2965010205431240,	0.0155951400782590},
        {0.0287564059530710,	0.2965010205431240,	0.1335581607035680,	0.5411844128002370,	0.0155951400782590}
    };
}
}
