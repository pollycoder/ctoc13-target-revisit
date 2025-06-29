#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include<math.h>
#include<vector>
#include<numeric>





//CTOC11相关

const int CityNum = 20;

const double g0 = 9.80665;//m/s^2
const double JD2S = 86400.0;
const double EarthRefRadius = 6378137.0;//m
const double Re             = 6378137.0;//m
const double Re_km          = 6378.137;//km
const double muEarth = 398600.4415e9;//m^3/s^2
const double omegaE     = 7.2921151467e-5;//rad/s
const double OmegaEarth = 7.2921151467e-5;//7.292158553e-5;//rad/s
const double mu_km_s = 398600.4415;//km^3/s^2
const double mu      = 398600.4415e9;//m^3/s^2
const double J2 =   1.08262668e-3; //J2摄动结果
const double Isp = 300.0;//s

const double alpha_G0 = 3.2310939479887431; //rad

// 单位为km
const double MU_EARTH = 398600.4415;  // km3/s2
const double we = 7.2921151467e-5;    // rad/s
const double Cd = 2.2;
const double s = 8e-6;                // km^2


//PI
const double DPI=3.1415926535897932384626433832795;
const double D2PI=6.283185307179586476925286766559;
const double D2R=0.017453292519943295769236907684886;
const double R2D=57.295779513082320876798154814105;
//Time
//const double JD2S=86400.0;
const double Y2JD=365.25;
//Constant
//const double mu=398600.4415e9;//m^3/s^2
//const double muEarth=398600.4415e9;//m^3/s^2
const double AU=1.49597870691e11;//m
//const double g0=9.80665;//m/s^2
const double EPSILON=1.0e-14;
const double ee=6.69437999014E-3;//地球椭圆形子午圈的偏心率平方
const double eep=6.73949674228E-3;//地球椭圆形子午圈的偏心率平方
//const double EarthRefRadius=6378137.0;//m
const double EarthPoleRadius=6356752.3142;//m
const double EarthStaRadius=42164172.367;//地球静止轨道半径
//const double OmegaEarth= 7.292158553e-5;//7.292158553e-5;//rad/s
//Unit
const double LUnit=AU;//一个天文单位
const double TUnit=sqrt(LUnit*LUnit*LUnit/mu);//半长轴为长度单位时，轨道周期在此时间单位下为2PI
//const double TUnit=365.25*86400;;//半长轴为长度单位时，轨道周期在此时间单位下为2PI
const double MuUnit=LUnit/TUnit*LUnit/TUnit*LUnit;
const double VUnit=LUnit/TUnit;
const double AUnit=VUnit/TUnit;
//Pro


const double MJDA0=2223;
const double MJDE0=3023;
const double SecondA0=192067200.00000000;
const double SecondE0=261187200.00000000;
const double JD2Second=86400.0;
const double mhNU=(mu/LUnit)*(TUnit/LUnit)*(TUnit/LUnit);
const double Second2Radian = 4.848136811095359935899141E-6;

//const double J2=1.0826269e-3;
//const double Re=6378137.0;//m
const double MJDJ2000 = 51544.5;
//const double mu=398600.4418e9;

const double JD_MJD  = 2400000.5;
const double TT_TAI = 32.184;
const double JC2JD= 36525.0;

//const double alpha_G0 = 1.747455428309031; //rad


//轨道机动相关
//const double Re_km = 6378.137;//km
//const double omegaE = 7.2921158553e-5;//rad/s
//const double mu_km_s = 398600.4415;//km^3/s^2
//const double Isp = 300.0;//s
//const double gamma = 2.5 / 180.0 * DPI;//rad
//const double alpha_G0 = 0.681733;//rad,论文场景



//CTOC13 相关
const double a_max = Re_km + 980.0;														    // Maximum SMA (km)
const double a0_min = Re_km + 520.0;														// Minimum SMA (km)
const double a_min = Re_km + 220.0;
const double e_scale = 0.003;																	// Maximum ECC

const double i_min = DPI / 2.0;
const double i_max = DPI;

const double dv_max = 1.0;

const double penalty = 1.0e6;

const int TreeNum = 8;
const int TargetNum = 21;


const int impNum = 8;


const double hmax = 1000.0;
const double hmin = 200.0;
const double imp_max = 2.0;

extern double t_anchor;
extern double t_gap;

// 第一个vector是脉冲星的编号，第二个vector是脉冲的次数
const std::vector<int> imp_sat = { 4, 5, 6, 7 };
const std::vector<int> imp_num = { 1, 1, 1, 1 };
extern std::vector<std::vector<std::vector<double>>> fixed_imp;

// 未完成的地面目标
const std::vector<int> ground_target_uncompleted = { 5, 16, 18 };




/*--------------------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------- 构型 -------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
// 数据库选的构型，8颗星全覆盖地面，80分
//const double sats_coe0[TreeNum][6] = { {7348.13700000001,	0.00100000000000042,	2.42600766,			3.64424747799999,	3.14159267466321,	0.502654824595922},
//									   {7298.13699999995,	0.00299999999999847,	2.42600766,			0.251327412299999,	1.5707963270021,	4.7752208329979},
//									   {7348.13699999998,	0.00100000000000203,	2.37364778300001,	1.193805208,		3.14159265358979,	2.38761041699917},
//									   {7298.13700000015,	0.00200000000000135,	2.164208272,		4.83805268699999,	6.28318528610617,	4.77522083300928},
//									   {7298.13700000008,	0.00199999999999933,	2.19911485800001,	6.094689748,		4.71238898000454,	2.13628300399546},
//									   {7300.13700000001,	0.00300000000000327,	2.129301687,		1.69646003300001,	0,					4.39822971500082},
//									   {7300.13700000001,	0.00300000000000187,	2.164208272,		3.26725636,			3.14159265358979,	4.5238934210007},
//									   {7300.1369999999,	0.00299999999998579,	2.37364778300001,	2.19911485799999,	3.14159265358979,	0.376991118404456} };

// 师弟们的构型，适用于大机动优化
//const double sats_coe0[TreeNum][6] = { {7220.860449638,		0.001175,				55.450306 * D2R,	153.159169 *  D2R,	68.298991 * D2R,	311.933840 * D2R},
//									   {7222.239050483,		0.000916,				55.280777 * D2R,	244.377245 *  D2R,	136.696849 * D2R,	48.486990 * D2R},
//									   {7218.206243053,		0.000400,				55.703358 * D2R,	62.185499 * D2R,	92.867924 * D2R,	122.298897 * D2R},
//									   {7215.001692415,		0.000867,				55.417036 * D2R,	331.082034 *  D2R,	111.759707 * D2R,	298.278435 * D2R},
//									   {7195.33169481176,	0.00137895034702387,	0.466644828766983,	3.57661573953852,	4.56945460441277,	6.11281774342053},
//									   {7196.94231594548,	0.000539730052198469,	0.465347943711639,	5.19980893281186,	0.0171067532199021,	0.732049838540783},
//									   {7198.47684002579,	0.0014840234629908,		0.470380281628898,	0.532534709005725,	3.84901229300942,	5.81928708190978},
//									   {7198.6102902429,	0.00133936973855959,	0.466335735561036,	2.06193691000045,	5.67963403118441,	0.543603078466579} };

// 大偏心率8颗星构型，机动优化已经复现出满分，脉冲3404m/s
const double sats_coe0[TreeNum][6] = { {7220.1873083705,	0.00014555827246772,	0.96959480984311,	2.6719021905352,	1.7372401063588,	4.8470637685766},
									   {7220.2529713729,	0.0040254132529364,		0.96481624891751,	4.2647308766549,	3.4295823654489,	6.000933555041},
									   {7218.8976541773,	0.0029084000693427,		0.97039802345484,	1.0843307086747,	3.1443601858976,	0.61471517678504},
									   {7216.3668294728,	0.001682757645501,		0.96758242339894,	5.7742711899074,	1.9914102044307,	5.1852416518209},
									   {7177.3170197352,	0.015560361518168,		0.49224353439439,	3.4553343550809,	5.6424963262825,	6.2446391554937},
									   {7108.8486429526,	0.028457351182438,		0.52189636807905,	5.1133971278817,	4.1555030830955,	2.663944787587},
									   {7274.2046679875,	0.011275823630897,		0.53196972751517,	0.55371780137885,	3.7246944101565,	6.1318317850596},
									   {7122.3513650047,	0.019440910345076,		0.45617922380697,	2.0101558037918,	3.97308136286,		4.0058692000737} };

//const double sats_coe0[TreeNum][6] = { {7220.1873083705,	0.00014555827246772,	0.96959480984311,	2.6719021905352,	1.7372401063588,	4.8470637685766},};




// 无机动12颗星满分构型
//const double sats_coe0[TreeNum][6] = { {7348.13700000001,	0.00100000000000042,	2.42600766,			3.64424747799999,	3.14159267466321,	0.502654824595922},
//									   {7298.13699999995,	0.00299999999999847,	2.42600766,			0.251327412299999,	1.5707963270021,	4.7752208329979},
//									   {7348.13699999998,	0.00100000000000203,	2.37364778300001,	1.193805208,		3.14159265358979,	2.38761041699917},
//									   {7298.13700000015,	0.00200000000000135,	2.164208272,		4.83805268699999,	6.28318528610617,	4.77522083300928},
//									   {7298.13700000008,	0.00199999999999933,	2.19911485800001,	6.094689748,		4.71238898000454,	2.13628300399546},
//									   {7300.13700000001,	0.00300000000000327,	2.129301687,		1.69646003300001,	0,					4.39822971500082},
//									   {7300.13700000001,	0.00300000000000187,	2.164208272,		3.26725636,			3.14159265358979,	4.5238934210007},
//									   {7300.1369999999,	0.00299999999998579,	2.37364778300001,	2.19911485799999,	3.14159265358979,	0.376991118404456},
//									   {7375.34264072713,	0.0,					2.9845130209103,	0.0317513281018904,	0.0,				6.21969192761202},
//									   {7375.34264072725,	0.0,					2.9845130209103,	3.60888342458644,	0.0,				0.191495396032018},
//									   {7375.00124841908,	0.0,					2.96705972839036,	5.89492854103316,	0.0,				0.114566563065866},
//									   {7375.00124841903,	0.0,					2.96705972839036,	0.944873606370856,	0.0,				0.0571920411461987} };

// 4颗星全机动覆盖地面
//const double sats_coe0[TreeNum][6] = { {7220.1873083705,	0.00014555827246772,	0.96959480984311,	2.6719021905352,	1.7372401063588,	4.8470637685766},
//								 {7220.2529713729,	0.0040254132529364,		0.96481624891751,	4.2647308766549,	3.4295823654489,	6.000933555041},
//								 {7218.8976541773,	0.0029084000693427,		0.97039802345484,	1.0843307086747,	3.1443601858976,	0.61471517678504},
//								 {7216.3668294728,	0.001682757645501,		0.96758242339894,	5.7742711899074,	1.9914102044307,	5.1852416518209} };


#endif