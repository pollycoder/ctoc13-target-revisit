#include "single_impluse.h"
extern int solve_cubic(double a, double b, double c, double d, double* xout);


//单脉冲机动(指定天数的最小脉冲解)
//输入：
//		m0:			初始质量，kg
//		t0:			初始（脉冲）时刻，s
//		rv0[6]:		卫星初始位置速度，km，km/s
//		lambda:		目标点的地面经度，rad
//		phi:		目标点的地面纬度，rad
//		Day:		设定的天数，计算范围(Day-1)~Day
//		branch:		升交段0或降交段1
//		sign:		需要的dv的符号：1正，-1负
//输出：
//		flag:		计算成功返回1，否则返回0
//		mf:			机动后卫星质量，kg
//		tf:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
				int& flag0, double& mf, double& tf, double* dv, int& NR,const int branch, const int sign)
{
	double lambda, phi;//目标点的经纬度，转化成弧度
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6],ocoe0[6];//初始轨道根数
	double c1, c2, k1, k2, k3, A, B, C, D;//中间参数
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//代表相关参数在J2下对时间的一阶导数
	double mcoe[6], ocoe[6];//最终结果的平瞬根转换
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//考虑J2卫星实际位置速度

	flag0 = 0;

	//rv2coe(flag, ocoe0, rv0, mu_km_s);//输入参数转化成轨道根数
	//ocoe0[0] *= 1000.0;
	//O2M(ocoe0, coe0, J2);
	//coe0[0] /= 1000.0;
	rv2coe(flag, coe0, rv0, mu_km_s);//输入参数转化成轨道根数

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);
	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	//lambda = lambda0 * D2R;//目标点的经纬度，转化成弧度
	//phi = phi0 * D2R;
	lambda = lambda0;
	phi = phi0;

	Nmax = floor(double(Day) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//最多转移圈数
	Nmin = floor((double(Day) - 1.0) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//最少转移圈数

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//初始时刻GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//注意范围[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//升交和降交两支
		int bra = branch;
		{
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//代表omega+ft，是最终轨道的参数，但是轨道面没变

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//求解关于(根号as)的三次方程，保留fabs(as-a0)最小的正根
			A = fs + 2.0 * DPI * double(NRi) - 2.0 * e0 * sin(fs) - 2.0 * k2 * sin(fs) - k1 + k3 + 2.0 * e0 * k1 * cos(fs) - M00;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			as = 0.0;
			delta_a = 1.0e10;
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					as = root[i] * root[i];
				}
			}
			//if (as < Re_km || as > Re_km + 1000.0)continue;

			//线性J2模型相关
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);

			//线性J2模型下的三次方程求解am，即平半长轴
			A = A - (M_J2_ + omega_J2_) * tG_J2;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG_J2 * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			am = 0.0;
			delta_a = 1.0e10;
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					am = root[i] * root[i];
				}
			}
			if (am < Re_km + 200.0)continue;

			//此时am已求出，且i,Omega不变，需计算其他平轨道根数（近似解）
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			//if (em < 0.0 || em > 1.0)continue;
			if (am * (1.0 - em) < Re_km + 200.0)continue;
			//if (am * (1.0 + em) > Re_km + 1000.0)continue;

			//平根转化为瞬根，由于平瞬根转换程序中使用的Re为m，在这里也临时改成m来计算
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			//计算初始时刻的切向脉冲
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0]/1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

			if ((double)sign * dv_mod < 0.0)continue;//新增脉冲方向的选择

			if (fabs(dv_min) > fabs(dv_mod))
			{
				dv_min = dv_mod;

				//C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(am, -3.5);
				//Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - em * em) / (1.0 - em * em);
				//tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);//精确的J2飞行时间
				tmin = tG_J2;
				NR = NRi;
				//branch = bra;
			}
		}
	}

	v0_mod = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
	for (int i = 0; i < 3; i++)iv[i] = rv0[i + 3] / v0_mod;
	for (int i = 0; i < 3; i++)dv[i] = dv_min * iv[i];

	tf = tmin;
	mf = m0 * exp(-fabs(dv_min * 1000.0) / g0 / Isp);

	if (mf < 300.0)
		flag0 = 0;
	else
		flag0 = 1;
		

};


//单脉冲机动(指定天数的最小脉冲解)
//输入：
//		m0:			初始质量，kg
//		t0:			初始（脉冲）时刻，s
//		rv0[6]:		卫星初始位置速度，km，km/s
//		lambda:		目标点的地面经度，rad
//		phi:		目标点的地面纬度，rad
//		Day:		设定的天数，计算范围(Day-1)~Day
//输出：
//		flag:		计算成功返回1，否则返回0
//		mf:			机动后卫星质量，kg
//		tf:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR) 
{
	double lambda, phi;//目标点的经纬度，转化成弧度
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6], ocoe0[6];//初始轨道根数
	double c1, c2, k1, k2, k3, A, B, C, D;//中间参数
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//代表相关参数在J2下对时间的一阶导数
	double mcoe[6], ocoe[6];//最终结果的平瞬根转换
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//考虑J2卫星实际位置速度

	flag0 = 0;

	//rv2coe(flag, ocoe0, rv0, mu_km_s);//输入参数转化成轨道根数
	//ocoe0[0] *= 1000.0;
	//O2M(ocoe0, coe0, J2);
	//coe0[0] /= 1000.0;
	rv2coe(flag, coe0, rv0, mu_km_s);//输入参数转化成轨道根数

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);
	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	//lambda = lambda0 * D2R;//目标点的经纬度，转化成弧度
	//phi = phi0 * D2R;
	lambda = lambda0;
	phi = phi0;

	Nmax = floor(double(Day) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//最多转移圈数
	Nmin = floor((double(Day) - 1.0) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//最少转移圈数

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//初始时刻GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//注意范围[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		for (int bra = 0; bra < 2; bra++)//升交和降交两支
		//int bra = branch;
		{
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//代表omega+ft，是最终轨道的参数，但是轨道面没变

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//求解关于(根号as)的三次方程，保留fabs(as-a0)最小的正根
			A = fs + 2.0 * DPI * double(NRi) - 2.0 * e0 * sin(fs) - 2.0 * k2 * sin(fs) - k1 + k3 + 2.0 * e0 * k1 * cos(fs) - M00;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			as = 0.0;
			delta_a = 1.0e10;
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					as = root[i] * root[i];
				}
			}
			if (as < Re_km || as > Re_km + 1500.0)continue;

			//线性J2模型相关
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);

			//线性J2模型下的三次方程求解am，即平半长轴
			A = A - (M_J2_ + omega_J2_) * tG_J2;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG_J2 * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			am = 0.0;
			delta_a = 1.0e10;
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					am = root[i] * root[i];
				}
			}
			if (am < Re_km + 200.0)continue;

			//此时am已求出，且i,Omega不变，需计算其他平轨道根数（近似解）
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			//if (em < 0.0 || em > 1.0)continue;
			if (am * (1.0 - em) < Re_km + 200.0)continue;
			if (am * (1.0 + em) > Re_km + 1500.0)continue;

			//平根转化为瞬根，由于平瞬根转换程序中使用的Re为m，在这里也临时改成m来计算
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			//计算初始时刻的切向脉冲
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

			//if ((double)sign * dv_mod < 0.0)continue;//新增脉冲方向的选择

			if (fabs(dv_min) > fabs(dv_mod))
			{
				dv_min = dv_mod;

				//C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(am, -3.5);
				//Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - em * em) / (1.0 - em * em);
				//tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);//精确的J2飞行时间
				tmin = tG_J2;
				NR = NRi;
				//branch = bra;
			}
		}
	}

	v0_mod = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
	for (int i = 0; i < 3; i++)iv[i] = rv0[i + 3] / v0_mod;
	for (int i = 0; i < 3; i++)dv[i] = dv_min * iv[i];

	tf = tmin;
	mf = m0 * exp(-fabs(dv_min * 1000.0) / g0 / Isp);

	//if (mf < 300.0)
	//	flag0 = 0;
	//else
	//	flag0 = 1;
	flag0 = 1;
}


void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const double dt,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, const int sign) {
	double lambda, phi;//目标点的经纬度，转化成弧度
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6], ocoe0[6];//初始轨道根数
	double c1, c2, k1, k2, k3, A, B, C, D;//中间参数
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//代表相关参数在J2下对时间的一阶导数
	double mcoe[6], ocoe[6];//最终结果的平瞬根转换
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//考虑J2卫星实际位置速度

	flag0 = 0;

	//rv2coe(flag, ocoe0, rv0, mu_km_s);//输入参数转化成轨道根数
	//ocoe0[0] *= 1000.0;
	//O2M(ocoe0, coe0, J2);
	//coe0[0] /= 1000.0;
	rv2coe(flag, coe0, rv0, mu_km_s);//输入参数转化成轨道根数

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);
	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	//lambda = lambda0 * D2R;//目标点的经纬度，转化成弧度
	//phi = phi0 * D2R;
	lambda = lambda0;
	phi = phi0;

	Nmax = floor(dt / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s))) + 1;//最多转移圈数
	Nmin = floor(dt / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//最少转移圈数

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//初始时刻GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//注意范围[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//升交和降交两支
		int bra = branch;
		{
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//代表omega+ft，是最终轨道的参数，但是轨道面没变

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_)) / omegaE;

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//求解关于(根号as)的三次方程，保留fabs(as-a0)最小的正根
			A = fs + 2.0 * DPI * double(NRi) - 2.0 * e0 * sin(fs) - 2.0 * k2 * sin(fs) - k1 + k3 + 2.0 * e0 * k1 * cos(fs) - M00;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			as = 0.0;
			delta_a = 1.0e10;
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					as = root[i] * root[i];
				}
			}
			//if (as < Re_km || as > Re_km + 1000.0)continue;

			//线性J2模型相关
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_) / (omegaE - Omega_J2_);

			//线性J2模型下的三次方程求解am，即平半长轴
			A = A - (M_J2_ + omega_J2_) * tG_J2;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG_J2 * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			am = 0.0;
			delta_a = 1.0e10;
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					am = root[i] * root[i];
				}
			}
			//if (am < Re_km + 200.0)continue;

			//此时am已求出，且i,Omega不变，需计算其他平轨道根数（近似解）
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			//if (em < 0.0 || em > 1.0)continue;
			//if (am * (1.0 - em) < Re_km + 200.0)continue;
			//if (am * (1.0 + em) > Re_km + 1000.0)continue;

			//平根转化为瞬根，由于平瞬根转换程序中使用的Re为m，在这里也临时改成m来计算
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			//计算初始时刻的切向脉冲
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

			if ((double)sign * dv_mod < 0.0)continue;//新增脉冲方向的选择

			if (fabs(dv_min) > fabs(dv_mod))
			//if(fabs(tG_J2 - dt < 0.7 * 3600.0))
			{
				dv_min = dv_mod;

				//C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(am, -3.5);
				//Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - em * em) / (1.0 - em * em);
				//tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);//精确的J2飞行时间
				tmin = tG_J2;
				NR = NRi;
				flag0 = 1;
				//branch = bra;
			}
		}
	}

	v0_mod = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
	for (int i = 0; i < 3; i++)iv[i] = rv0[i + 3] / v0_mod;
	for (int i = 0; i < 3; i++)dv[i] = dv_min * iv[i];

	tf = tmin;
	//mf = m0 * exp(-fabs(dv_min * 1000.0) / g0 / Isp);

	/*if (mf < 300.0)
		flag0 = 0;
	else*/
}


double obj_func(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {	
	// Parameters
	double* para = static_cast<double*> (f_data);
	const double time_stamp = para[10];
	
	// Initialize the variables
	double dv[3] = { 0.0, 0.0, 0.0 };
	double rvf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double tf = 0.0;

	get_score_data(X, para, dv, rvf, tf);
	return V_Norm2(dv, 3);
}



void get_score_data(const std::vector<double>& X, const double* para, double* dv, double* rvf, double& tf) {
	const double hour = 3600.0;
	const double t0 = para[0];
	tf = para[1];
	const double rv0[6] = { para[2], para[3], para[4], para[5], para[6], para[7] };
	const double lambda0 = para[8];
	const double phi0 = para[9];
	const double time_stamp = para[10];
	const double dv0[3] = { para[11], para[12], para[13] };

	if(tf < 0) {
		//std::cout << "tf = " << tf << std::endl;
		for (int j = 0; j < 3; j++) {
			dv[j] = penalty;
			rvf[j] = penalty; rvf[j + 3] = penalty;
			tf = penalty;
		}
		return;
	}

	// Perturbation
	dv[0] = dv0[0] + (X[0] - 0.5) * 1e3;
	dv[1] = dv0[1] + (X[1] - 0.5) * 1e3;
	dv[2] = dv0[2] + (X[2] - 0.5) * 1e3;
	tf += (X[3] - 0.5) * 2 * hour;															// Time: ±0.5h

	// Observing time: if the time is too far away from the time stamp, return penalty
	/*if (fabs(time_stamp - tf) > 0.7 * hour) {
		for (int j = 0; j < 3; j++) {
			dv[j] = penalty;
			rvf[j] = penalty; rvf[j + 3] = penalty;
			tf = penalty;
		}
		return;
	}*/

	// Add impulse and propagate
	double v0[3] = { rv0[3], rv0[4], rv0[5] };
	double v0_norm = V_Norm2(v0, 3);
	V_Add(v0, v0, dv, 3);
	double rv1[6] = { rv0[0], rv0[1], rv0[2], v0[0], v0[1], v0[2] };
	propagate_j2(rv1, rvf, t0, tf);

	// Visibility: if the target is not visible, return penalty
	double R_target[3];
	double geo[2] = { phi0, lambda0 };
	Geodetic2J2000(geo, R_target, tf);
	double half_cone_angle = 19.5 * D2R;
	bool ifVisible = is_target_visible(rvf, R_target, half_cone_angle);
	if (!ifVisible) { 
		for (int j = 0; j < 3; j++) {
			dv[j] = penalty;
			rvf[j] = penalty; rvf[j + 3] = penalty;
			tf = penalty;
		}
		return; 
	}
	/*else {
		std::cout << "Succeed." << std::endl;
		for (int j = 0; j < 6; j++) {
			std::cout << "rvf[" << j << "] = " << rvf[j] << std::endl;
		}
	}*/
}


void shooting_target2target(const double t0, const double* rv0, const double time_stamp, const double lambda0, const double phi0, double& tf, double* dv, double* rvf, const int& branch) {
	const double m0 = 1.0e6;																				// Useless, randomly assigned one valid value
	double mf = 0.0; int flag = -1; int NR = -1;

	// Initial value: single_imp
	// All the 4 situations should be taken into account: 
	// +dv, ascend; -dv, ascend; +dv, descend; -dv, descend
	// Choose the one without penalty (mostly just one valid solution)
	// tf: final time, not flight time, therefore t0 should be added after single_impulse provides an initial value
	double impulse_temp = 1.0e6;
	double impulse = 1.0e6;
	double dt = time_stamp - t0;
	for (int dv_sign = -1; dv_sign < 2; dv_sign += 2) {
		single_imp(m0, t0, rv0, lambda0, phi0, dt, flag, mf, tf, dv, NR, branch, dv_sign);
		if (flag == 0) {
			tf = penalty;
			for (int i = 0; i < 3; i++) {
				dv[i] = penalty;
				rvf[i] = penalty;
				rvf[i + 3] = penalty;
			}
			continue;
		}
		tf += t0;																	
		double f_data[14] = { t0, tf, rv0[0], rv0[1], rv0[2], rv0[3], rv0[4], rv0[5], lambda0, phi0, time_stamp, dv[0], dv[1], dv[2] };
		std::vector<double> X = { 0.5, 0.5, 0.5, 0.5 };
		nlopt_main(obj_func, f_data, X, impulse_temp, 4);

		if (impulse_temp < impulse) {
			get_score_data(X, f_data, dv, rvf, tf);
			impulse = impulse_temp;
			/*std::cout << "tf = " << tf << std::endl;
			std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
			std::cout << "NR = " << NR << std::endl;*/
			return;
		}
	}
}


//CTOC13：利用J2Lambert问题，进行单次脉冲修正
//与张刚论文作用类似，但是脉冲并不是近似切向，在一个轨道高度范围内近似选择脉冲最低且能观测到地面目标的
//输入：
//		RV0[6]：初始位置速度
//		t0：初始时刻
//		tf：终端时刻
//输出：
//		dv[3]：最终的脉冲
//		RVf[6]：终端位置速度
//		flag：求解成功返回1，求解失败返回0    
//		h：最终采用的高度   
void single_imp(double* dv, double* RVf, int& flag, const double* RV0, const double& t0, const double& tf, const double& lambda0, const double& phi0, const double& h) {
	double R0[3] = { RV0[0], RV0[1], RV0[2] };
	double V0[3] = { RV0[3], RV0[4], RV0[5] };
	double geogetic_Target[2] = { phi0, lambda0 };
	double R_Target[3];
	Geodetic2J2000(geogetic_Target, R_Target, tf);

	int flag_temp; flag = 0;
	
	double dv_norm = 1.0e10;
	double dv_temp = 1.0e10;
	double v1temp[3], dvtemp[3], RVf_temp[6], RV1_temp[6];
	
	double R = Re_km + h;
	double Rf_temp[3];
	V_Multi(Rf_temp, R_Target, R / Re_km, 3);
	for (int i = 0; i < 3; i++) RVf_temp[i] = Rf_temp[i];
	J2Lambert_short(flag_temp, RV1_temp, RVf_temp, RV0, tf - t0, mu_km_s);

	bool ifvisible = is_target_visible(RVf_temp, R_Target, 20.0 * D2R);
	if (!ifvisible) flag_temp = 0;

	//如果不成功则惩罚
	if (flag_temp != 1) {
		//h = penalty;
		for (int j = 0; j < 3; j++) {
			RVf[j] = penalty;
			RVf[j + 3] = penalty;
			dv[j] = penalty;
		}

		flag = 0;
		return;
	}

	//如果成功则计算脉冲
	for (int i = 0; i < 3; i++) v1temp[i] = RV1_temp[i + 3];
	V_Minus(dvtemp, v1temp, V0, 3);
	dv_temp = V_Norm2(dvtemp, 3);
	if (dv_temp < dv_norm) {
		memcpy(RVf, RVf_temp, 6 * sizeof(double));
		memcpy(dv, dvtemp, 3 * sizeof(double));

		flag = 1;
		dv_norm = dv_temp;
	}
}

//CTOC13：利用J2Lambert，优化覆盖目标点的轨道高度
//优化变量（被摄动项）：
//		h：轨道高度
//		(lambda, phi)：实际经过的地面经纬度，先扰动时间，再进一步扰动经纬度，扰动量180°
//		tf：实际的终末时刻，从3h开始扰动，区间为[0,6h]
//优化指标：
//		dv：脉冲大小
//参数：
//		t0：初始时刻
//		RV0[6]：初始位置速度
//		tf：优化前的终末时刻
//		h0：优化前的轨道高度
//		target_id：要观测的目标序号（0-20）
//约束：
//		保证目标可见
//		轨道高度200-1000km
//格式与PSO对应的Obj_func完全一致，我们只需要获取最终的优化变量X
//不同于之前每次都要先写get_value再包装，这是为了childnode函数看起来思路更清晰
void perturbation(double& h, double& tf, double& lambda, double& phi, const std::vector<double>& X, const double& h0, const int& id) {
	const double double_h_pert = 100.0;
	const double double_angle_pert = 12.0 * D2R;
	h = h0 + (X[0] - 0.5) * double_h_pert;
	tf += (X[3] - 0.5) * 43200.0;								// 扰动经纬度之前，先扰动时间

	// 获取时间扰动后的经纬度（对地面目标无影响，主要影响海上目标）
	double target_geo[2];
	get_target_geogetic(id, tf, target_geo);
	const double lambda0 = target_geo[1];
	const double phi0 = target_geo[0];

	// 对时间扰动后的经纬度再做角度扰动（对所有目标都有影响）
	lambda = lambda0 + (X[1] - 0.5) * double_angle_pert;
	phi = phi0 + (X[2] - 0.5) * double_angle_pert;
}


double obj_func_shooting(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	double* para = static_cast<double*>(f_data);

	// 初值
	const double t0 = para[0];
	double tf = para[1];
	const double h0 = para[2];
	const double RV0[6] = { para[3], para[4], para[5], para[6], para[7], para[8] };
	const int id = static_cast<int>(para[9]);

	double lambda, phi, h;
	perturbation(h, tf, lambda, phi, X, h0, id);

	if (tf < t0) {
		//std::cout << "tf < t0，不符合要求" << std::endl;
		//std::cout << std::endl;
		return penalty;
	}

	double dv[3] = { 0.0, 0.0, 0.0 };
	double impulse = penalty;
	double RVf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	int flag = -1;
	single_imp(dv, RVf, flag, RV0, t0, tf, lambda, phi, h);

	if(h < 220.0 || h > 980.0) { 
		//std::cout << "高度不符合要求, h = " << h << " km" << std::endl;
		//std::cout << std::endl;
		return penalty; 
	}
	if(flag != 1) { 
		//std::cout << "J2 Lambert 求解失败" << std::endl;
		//std::cout << std::endl;
		return penalty; 
	}

	double target_R[3];
	get_target_R(id, tf, target_R);
	bool ifVisible = is_target_visible(RVf, target_R, 19.5 * D2R);
	if(!ifVisible) { 
		//std::cout << "打靶后目标不可见" << std::endl;
		//std::cout << std::endl;
		return penalty; 
	}

	impulse = V_Norm2(dv, 3);
	
	/*std::cout << "X = " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << std::endl;
	std::cout << "ID = " << id << std::endl;
	std::cout << "dv = " << impulse << std::endl;
	std::cout << std::endl;*/
	return impulse;
}


void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const double& h0, const int& target_id) {
	double f_data[10] = { t0, tf, h0, RV0[0], RV0[1], RV0[2], RV0[3], RV0[4], RV0[5], target_id };

	double impulse = 0.0;
	std::vector<double> X = { 0.5, 0.5, 0.5, 0.5 };
	nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 100);		//不输出

	double lambda, phi, h;
	perturbation(h, tf, lambda, phi, X, h0, target_id);

	single_imp(dv, RVf, flag, RV0, t0, tf, lambda, phi, h);
	if (flag != 1) {
		for (int j = 0; j < 3; j++) {
			dv[j] = penalty;
			RVf[j] = penalty; RVf[j + 3] = penalty;
		}
		tf = penalty;
	}
}