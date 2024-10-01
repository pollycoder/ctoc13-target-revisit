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
const double e_scale = 0.005;																	// Maximum ECC

const double i_min = DPI / 2.0;
const double i_max = DPI;

const double penalty = 1.0e6;

const int TreeNum = 1;
const int TargetNum = 21;

const int GapNum = 26;

#endif