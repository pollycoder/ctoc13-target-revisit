#ifndef CONSTANT
#define CONSTANT
#include<math.h>

const int TreeNum = 5;
const int CityNum = 20;

const double Day2Second = 86400.0;      //天转s
const double mu = 398600.4418e9;	    //地球引力常数(m^3/s^2)
const double req = 6378137.0;			//地球半径(m)
const double J2 = 1.08262668e-3;		//J2摄动
const double MJD_Init= 23467.0;         //初始时刻

const double alpha_G0 = 1.747455428309031; //rad
const double Isp = 300.0;//s


const double Re_km = 6378.137;//km
const double omegaE = 7.2921158553e-5;//rad/s
const double mu_km_s = 398600.4415;//km^3/s^2

//PI
const double pi = 3.1415926535897932384626433832795;	//圆周率
const double DPI = 3.1415926535897932384626433832795;
const double D2PI = 6.283185307179586476925286766559;
const double D2R = 0.017453292519943295769236907684886;
const double R2D = 57.295779513082320876798154814105;

const double gamma = 2.5 / 180.0 * DPI;//rad

const double EPSILON = 1.0e-14;
const double ee = 6.69437999014E-3;//地球椭圆形子午圈的偏心率平方
const double eep = 6.73949674228E-3;//地球椭圆形子午圈的偏心率平方


//Time
const double JD2S = 86400.0;
const double Y2JD = 365.25;
//Constant

const double muEarth = 398600.4415e9;//m^3/s^2
const double AU = 1.49597870691e11;//m
const double g0 = 9.80665;//m/s^2



const double EarthRefRadius = 6378137.0;//m
const double EarthPoleRadius = 6356752.3142;//m
const double EarthStaRadius = 42164172.367;//地球静止轨道半径
const double OmegaEarth = 7.2921158553e-5;//7.292158553e-5;//rad/s
//Unit
const double LUnit = AU;//一个天文单位
const double TUnit = sqrt(LUnit * LUnit * LUnit / mu);//半长轴为长度单位时，轨道周期在此时间单位下为2PI
//const double TUnit=365.25*86400;;//半长轴为长度单位时，轨道周期在此时间单位下为2PI
const double MuUnit = LUnit / TUnit * LUnit / TUnit * LUnit;
const double VUnit = LUnit / TUnit;
const double AUnit = VUnit / TUnit;
//Pro


const double MJDA0 = 2223;
const double MJDE0 = 3023;
const double SecondA0 = 192067200.00000000;
const double SecondE0 = 261187200.00000000;
const double JD2Second = 86400.0;
const double mhNU = (mu / LUnit) * (TUnit / LUnit) * (TUnit / LUnit);
const double Second2Radian = 4.848136811095359935899141E-6;


const double Re = 6378137.0;//m
const double MJDJ2000 = 51544.5;
//const double mu=398600.4418e9;

const double JD_MJD = 2400000.5;
const double TT_TAI = 32.184;
const double JC2JD = 36525.0;




//轨道机动相关





#endif

