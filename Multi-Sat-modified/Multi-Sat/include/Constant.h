#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include<math.h>






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

const int TreeNum = 12;
const int TargetNum = 21;

const int impNum = 4;

// 师弟们的构型，适用于大机动优化
//const double sats_coe0[TreeNum][6] = { {7220.860449638,		0.001175,				0.96779,			2.67313,			1.19204,			5.44427},
//									   {7222.239050483,		0.000916,				0.964832,			4.26519,			2.38581,			0.846258},
//									   {7218.206243053,		0.000400,				0.972207,			1.08534,			1.62085,			2.13452},
//									   {7215.001692415,		0.000867,				0.96721,			5.77847,			1.95057,			5.20594},
//									   {7195.33169481176,	0.00137895034702387,	0.466644828766983,	3.57661573953852,	4.56945460441277,	6.11281774342053},
//									   {7196.94231594548,	0.000539730052198469,	0.465347943711639,	5.19980893281186,	0.0171067532199021,	0.732049838540783},
//									   {7198.47684002579,	0.0014840234629908,		0.470380281628898,	0.532534709005725,	3.84901229300942,	5.81928708190978},
//									   {7198.6102902429,	0.00133936973855959,	0.466335735561036,	2.06193691000045,	5.67963403118441,	0.543603078466579} };


// 暂时用无机动构型优化
const double sats_coe0[TreeNum][6] = { {7348.13700000001,	0.00100000000000042,	2.42600766,			3.64424747799999,	3.14159267466321,	0.502654824595922},
									   {7298.13699999995,	0.00299999999999847,	2.42600766,			0.251327412299999,	1.5707963270021,	4.7752208329979},
									   {7348.13699999998,	0.00100000000000203,	2.37364778300001,	1.193805208,		3.14159265358979,	2.38761041699917},
									   {7298.13700000015,	0.00200000000000135,	2.164208272,		4.83805268699999,	6.28318528610617,	4.77522083300928},
									   {7298.13700000008,	0.00199999999999933,	2.19911485800001,	6.094689748,		4.71238898000454,	2.13628300399546},
									   {7300.13700000001,	0.00300000000000327,	2.129301687,		1.69646003300001,	0,					4.39822971500082},
									   {7300.13700000001,	0.00300000000000187,	2.164208272,		3.26725636,			3.14159265358979,	4.5238934210007},
									   {7300.1369999999,	0.00299999999998579,	2.37364778300001,	2.19911485799999,	3.14159265358979,	0.376991118404456},
									   {7375.34264072713,	0.0,					2.9845130209103,	0.0317513281018904,	0.0,				6.21969192761202},
									   {7375.34264072725,	0.0,					2.9845130209103,	3.60888342458644,	0.0,				0.191495396032018},
									   {7375.00124841908,	0.0,					2.96705972839036,	5.89492854103316,	0.0,				0.114566563065866},
									   {7375.00124841903,	0.0,					2.96705972839036,	0.944873606370856,	0.0,				0.0571920411461987} };

#endif