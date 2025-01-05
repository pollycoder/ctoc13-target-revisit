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

	Nmax = floor(10800.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//最多转移圈数(3h)
	Nmin = 0;//floor((double(Day) - 1.0) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//最少转移圈数

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
			if (as < Re_km || as > Re_km + 1000.0)continue;

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
//TODO：去掉里面branch和NR的遍历，也不要求脉冲排序，求解成功即可扩展
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
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch) 
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

	rv2coe(flag0, coe0, rv0, mu_km_s);//输入参数转化成轨道根数

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

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//初始时刻GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//注意范围[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	//for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//升交和降交两支
		int bra = branch;
		NRi = NR;
		{
			// ZZ 这里错了
			if (phi > inc0 + 0.1) {
				flag0 = 0;
				return;
			}
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
			if (root_num < 1) {
				flag0 = 0;
				return;
			}
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					as = root[i] * root[i];
				}
			}

			/*if (as < Re_km) {
				flag0 = 0;
				return;
			}
			*/

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
			
			//此时am已求出，且i,Omega不变，需计算其他平轨道根数（近似解）
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);

			double peri_m = am * (1 - em);
			double apo_m = am * (1 + em);

			if (em < 0.0 || em > 1.0) {
				flag0 = 0;
				return;
			}
			/*if (peri_m  < Re_km + 180.0) {
				flag0 = 0;
				return;
			}
			if (apo_m > Re_km + 1100.0) {
				flag0 = 0;
				return;
			}*/
			if (fabs(am - 0.0) < 1.0e-5) {
				flag0 = 0;
				return;
			}

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
			double temp1 = mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0));
			double temp2 = mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00));
			if (temp1 < 0 || temp2 < 0) {
				flag0 = 0;
				return;
			}
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));
			dv_min = dv_mod;

			// 判断是否为nan
			if (dv_min != dv_min) {
				flag0 = 0;
				return;
			}

			tmin = tG_J2;
		}
	}

	v0_mod = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
	for (int i = 0; i < 3; i++)iv[i] = rv0[i + 3] / v0_mod;
	for (int i = 0; i < 3; i++)dv[i] = dv_min * iv[i];

	tf = tmin;
	mf = m0 * exp(-fabs(dv_min * 1000.0) / g0 / Isp);

	flag0 = 1;
}


void calculate_tg_2body(const int Day, double lambda, double phi, double inc0, double Omega0, double& tG, double alphaG0_, int bra)
{
	double c1;
	double c2;
	double alphaGt;
	double omega_ft;

	omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//代表omega+ft，是最终轨道的参数，但是轨道面没变

	c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
	c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
	alphaGt = Omega0 + atan2(c1, c2);
	while (alphaGt < alphaG0_)alphaGt += D2PI;
	while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
	tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;
}


//Tg 是经过多长时间，不是时刻
void modify_tg(const double t0, const double* rv0, double phi, double& tG)
{
	//再次修正时刻
	double rv_temp[6];
	memcpy(rv_temp, rv0, 6 * sizeof(double));
	propagate_linearJ2(rv0, rv_temp, t0, tG + t0);

	int flag_coe;
	double coe[6];
	rv2coe(flag_coe, coe, rv_temp, mu * 1e-9);  //求解轨道根数

	double u0 = coe[4] + coe[5]; u0 = NiceAngle(u0);//相角//范围在0 到2pi
	double u1 = asin(sin(phi) / sin(coe[2]));
	double inc_lower90 = coe[2]; if (inc_lower90 > DPI / 2.0)  inc_lower90 = DPI - inc_lower90;
	if (phi > inc_lower90 + 0.06)
	{
		return;
	}
	else if(phi > inc_lower90 - 0.06 )
	{
		u1 = asin(sin(phi) / sin(phi )); //inc_lower90 + 0.00 
	}
	else
	{
		u1 = asin(sin(phi) / sin(coe[2]));
	}


	
	double u2 = DPI - u1;
	//范围不能相差pi 
	double du1 = NiceAngle(u1 - u0);
	double du2 = NiceAngle(u2 - u0);
	if (du1 > DPI) du1 = du1 - 2 * DPI;
	if (du2 > DPI) du2 = du2 - 2 * DPI;

	//start时刻的位置速度
	double r_start_staellite[3] = { rv_temp[0],rv_temp[1],rv_temp[2] };
	double v_start_staellite[3] = { rv_temp[3],rv_temp[4],rv_temp[5] };

	//寻找星下点纬度第一次与目标重合的两个时间
	double _omega_ = V_Norm2(v_start_staellite, 3) / V_Norm2(r_start_staellite, 3);
	double t1 = du1 / _omega_;
	double t2 = du2 / _omega_;
	double t = fabs(t1) < fabs(t2) ? t1 : t2;
	double  t_access = tG + t;

	if (t_access > 0)
    	tG = t_access;
}

//单脉冲机动(指定天数的最小脉冲解:纯二体Lambert求解)
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
//		tG:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
void single_imp_zzmodified(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tG, double* dv, const int& branch)
{
	double lambda, phi;//目标点的经纬度，转化成弧度
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6], ocoe0[6];//初始轨道根数
	double c1, c2;//中间参数
	double alphaG0_, alphaGt;//GMST
	double omega_ft;
	int flag;

	flag0 = 0;

	rv2coe(flag0, coe0, rv0, mu_km_s);//输入参数转化成轨道根数

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);

	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	lambda = lambda0;
	phi = phi0;

	alphaG0_ = alpha_G0 + omegaE * t0;//初始时刻GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//注意范围[0, 2PI]


	int bra = branch;

	double inc_lower90 = inc0;
	if (inc_lower90 > DPI / 2.0) inc_lower90 = DPI - inc_lower90;

	if(phi > inc0 + 0.06)
	{
		flag0 = 0;
		return;
	}

	//一半张角来计算
	if (phi > inc0 -0.06) {
		//ZZ TODO: 需要考虑86400 跳天的情况

		double tG_temp;
		calculate_tg_2body(Day, lambda, phi, inc_lower90+0.06, Omega0, tG_temp, alphaG0_, bra);

		if (86400 - tG_temp < 7200)  tG_temp = tG_temp - 86400.0;

		double tG_temp2;
		calculate_tg_2body(Day, lambda, phi, phi, Omega0, tG_temp2, alphaG0_, bra);
		if (86400 - tG_temp2 < 7200)  tG_temp2 = tG_temp2 - 86400.0;

		tG = (tG_temp + tG_temp2) / 2.0;
		if (tG < 0.0) tG = std::max(tG_temp, tG_temp2) / 1.2;
	}
	else 
	{
		double tG_temp;
		calculate_tg_2body(Day, lambda, phi, inc_lower90 + 0.06, Omega0, tG_temp, alphaG0_, bra);
		if (86400 - tG_temp < 7200)  tG_temp = tG_temp - 86400.0;
		double tG_temp2;
		calculate_tg_2body(Day, lambda, phi, inc_lower90 - 0.06, Omega0, tG_temp2, alphaG0_, bra);
		if (86400 - tG_temp2 < 7200)  tG_temp2 = tG_temp2 - 86400.0;
		tG = (tG_temp + tG_temp2) / 2.0;
		if (tG < 0.0) tG = std::max(tG_temp, tG_temp2);
	}

	modify_tg(t0, rv0, phi, tG);

	if(tG > 20000.0)
	{
		flag0 = 0;
		return;
	}

	 //计算精确过点时刻的R向量
	double Geodetic[2] = { phi0,lambda0 };
	double RV_Target[6];
	Geodetic2J2000(Geodetic, RV_Target, tG + t0);
	V_Multi(RV_Target, RV_Target, a0 / V_Norm2(RV_Target, 3), 3);

	//////////////////////////////////////计算Lambert问题////////////////////////
	double vt0[6], am, vt1[3], v0[3], v1[3], a, e, unith[3], dtemp0, dtemp1, dvmin = 1.0e10;
	int  Nmax;
	flag = 0;
	
	V_Cross(unith, rv0, &rv0[3]);
	dtemp0 = V_Norm2(unith, 3);
	if (dtemp0 > 0.0) for (int i = 0; i < 3; i++) unith[i] /= dtemp0;
	else
	{
		unith[0] = 0.0;
		unith[1] = 0.0;
		unith[2] = 1.0;
	}
	for (int i = 0; i < 3; i++) vt0[i] = rv0[i] - RV_Target[i];
	am = (V_Norm2(rv0, 3) + V_Norm2(RV_Target, 3) + V_Norm2(vt0, 3)) / 4.0; //最小椭圆转移轨道半长轴
	Nmax = (int)floor(tG / (6.283185307179586476925286766559 * sqrt(am * am * am / mu_km_s))) + 1; //最多可能的转移圈数
	int flag_total = 0;

	double dt_min = 0.0;

	for (int n = 0; n <= Nmax; n++) //圈数从0到最大可能数穷举，保留特征速度最小的结果
	{
		for (int j = 0; j < 2; j++) //对于多圈问题，有左右分枝两个解
		{
			double dt = 0.0;
			//for(double dt = -2000.0; dt < 2010.0 ; dt+=100.0)
			{
				
				//计算精确过点时刻的R向量
				//double Geodetic[2] = { phi0,lambda0 };
				//double RV_Target[6];
				if(tG + dt < 50.0) continue;

				double d_phi = -0.0;
				//for (double d_phi = - 0.07; d_phi < 0.071 ; d_phi+= 0.03)
				{
					Geodetic[0] =  phi0 + d_phi ;

					double d_lamda = 0.0;
					//for (double d_lamda = -0.07; d_lamda < 0.071; d_lamda += 0.03)
					{
						Geodetic[1] =  lambda0 + d_lamda ;

						//for(double R_norm = Re_km + 500.0; R_norm < Re_km + 1000.1; R_norm += 10.0)
						{
							Geodetic2J2000(Geodetic, RV_Target, tG + t0 + dt);
							V_Multi(RV_Target, RV_Target, a0 / V_Norm2(RV_Target, 3), 3);

							lambert(v0, v1, a, e, rv0, RV_Target, tG, unith, flag0, mu_km_s, 0, n, j, 60);

							if (flag0 != 1 || e >= 1.0) continue;//排除双曲和抛物轨道
							for (int i = 0; i < 3; i++)
							{
								vt0[i] = v0[i] - rv0[3 + i];
								vt1[i] = RV_Target[3 + i] - v1[i];
							}
							dtemp0 = sqrt(vt0[0] * vt0[0] + vt0[1] * vt0[1] + vt0[2] * vt0[2]); // V_Norm2(vt0, 3);
							dtemp1 = sqrt(vt1[0] * vt1[0] + vt1[1] * vt1[1] + vt1[2] * vt1[2]); // V_Norm2(vt1, 3);

							if (dtemp0 < dvmin)
							{
								dt_min = dt;
								flag_total = 1;
								flag = 1;
								dvmin = dtemp0;
								memcpy(dv, vt0, 3 * sizeof(double));
							}
						}

					}

				}

			}
		}
	}

	tG = tG + dt_min; //此处很关键，不能变成tG + dt_min，可以通过优化修正
	mf = m0 * exp(-fabs(dvmin * 1000.0) / g0 / Isp);
	if (flag_total == 1)
	{
		flag0 = 1;
	}
	else
	{
		flag0 = 0;
	}
}

//单脉冲机动(指定天数的最小脉冲解:纯二体Lambert求解)
void lambert_dv(int& flag, const double t0, const double* rv0,  double* dv,  double tf, int id)
{
	flag = 0;
	if (tf - t0 < 10.0) return;

	double coe0[6];
	int flag_orbit;
	rv2coe(flag_orbit, coe0, rv0, mu_km_s);//输入参数转化成轨道根数
	double a0 = coe0[0];

	double Geodetic[2];
	get_target_geogetic(id, tf, Geodetic);

	double RV_Target[6];
	Geodetic2J2000(Geodetic, RV_Target, tf);
	V_Multi(RV_Target, RV_Target, a0 / V_Norm2(RV_Target, 3), 3);

	double vt0[6], am, vt1[3], v0[3], v1[3], a, e, unith[3], dtemp0, dtemp1;
	double dvmin = 1.0e10;
	int  Nmax;

	V_Cross(unith, rv0, &rv0[3]);
	dtemp0 = V_Norm2(unith, 3);
	if (dtemp0 > 0.0) for (int i = 0; i < 3; i++) unith[i] /= dtemp0;
	else
	{
		unith[0] = 0.0;
		unith[1] = 0.0;
		unith[2] = 1.0;
	}
	for (int i = 0; i < 3; i++) vt0[i] = rv0[i] - RV_Target[i];
	am = a0; //最小椭圆转移轨道半长轴
	Nmax = (int)floor((tf-t0) / (6.283185307179586476925286766559 * sqrt(am * am * am / mu_km_s))) + 1; //最多可能的转移圈数


	for (int n = 0; n <= Nmax; n++) //圈数从0到最大可能数穷举，保留特征速度最小的结果
	{
		for (int j = 0; j < 2; j++) //对于多圈问题，有左右分枝两个解
		{
			
			for(int way = 0; way < 2; way++)
			{
				Geodetic2J2000(Geodetic, RV_Target, tf );
				V_Multi(RV_Target, RV_Target, a0 / V_Norm2(RV_Target, 3), 3);
				int flag0;
				lambert(v0, v1, a, e, rv0, RV_Target, tf  - t0, unith, flag0, mu_km_s, way, n, j, 200);

				if (flag0 != 1 || e >= 1.0) continue;//排除双曲和抛物轨道
				for (int i = 0; i < 3; i++)
				{
					vt0[i] = v0[i] - rv0[3 + i];
					vt1[i] = RV_Target[3 + i] - v1[i];
				}
				dtemp0 = sqrt(vt0[0] * vt0[0] + vt0[1] * vt0[1] + vt0[2] * vt0[2]); // V_Norm2(vt0, 3);
				dtemp1 = sqrt(vt1[0] * vt1[0] + vt1[1] * vt1[1] + vt1[2] * vt1[2]); // V_Norm2(vt1, 3);

				if (dtemp0 < dvmin)
				{
					flag = 1;
					dvmin = dtemp0;
					memcpy(dv, vt0, 3 * sizeof(double));
				}
			}

		}
	}
	
}


void single_imp_ship(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch)
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

	double lon1 = 53.558 * D2R;
	double lon2 = 77.963 * D2R;
	double omega_ship = (lon2 - lon1) / 2.0 / 86400.0;


	flag0 = 0;

	rv2coe(flag0, coe0, rv0, mu_km_s);//输入参数转化成轨道根数

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

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//初始时刻GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//注意范围[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	//for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//升交和降交两支
		int bra = branch;
		NRi = NR;
		{
			if (phi > inc0) {
				flag0 = 0;
				return;
			}
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//代表omega+ft，是最终轨道的参数，但是轨道面没变

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE + omega_ship);

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//求解关于(根号as)的三次方程，保留fabs(as-a0)最小的正根
			A = fs + 2.0 * DPI * double(NRi) - 2.0 * e0 * sin(fs) - 2.0 * k2 * sin(fs) - k1 + k3 + 2.0 * e0 * k1 * cos(fs) - M00;
			B = 0.0;
			C = (2.0 * k2 * sin(fs) + k1 - k3 - 2.0 * e0 * k1 * cos(fs)) * a0;
			D = -tG * sqrt(mu_km_s);
			root_num = solve_cubic(A, B, C, D, root);
			as = 0.0;
			delta_a = 1.0e10;
			if (root_num < 1) {
				flag0 = 0;
				return;
			}
			for (int i = 0; i < root_num; i++)
			{
				if (root[i] > 0 && delta_a > fabs(root[i] * root[i] - a0))
				{
					delta_a = fabs(root[i] * root[i] - a0);
					as = root[i] * root[i];
				}
			}

			/*if (as < Re_km) {
				flag0 = 0;
				return;
			}*/


			//线性J2模型相关
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE + omega_ship - Omega_J2_);

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

			//此时am已求出，且i,Omega不变，需计算其他平轨道根数（近似解）
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			if (em < 0.0 || em > 1.0) {
				flag0 = 0;
				return;
			}
			/*if (am < Re_km + 200.0) {
				flag0 = 0;
				return;
			}
			if (am > Re_km + 1000.0) {
				flag0 = 0;
				return;
			}*/

			//平根转化为瞬根，由于平瞬根转换程序中使用的Re为m，在这里也临时改成m来计算
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			/*if (ocoe[0] / 1000.0 * (1 - ocoe[1]) < Re_km + 200.0) {
				flag0 = 0;
				return;
			}

			if (ocoe[0] / 1000.0 * (1 + ocoe[1]) > Re_km + 1000.0) {
				flag0 = 0;
				return;
			}*/

			//计算初始时刻的切向脉冲
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));
			dv_min = dv_mod;

			// 判断是否为nan
			if (dv_min != dv_min) {
				flag0 = 0;
				return;
			}

			tmin = tG_J2;
		}
	}

	v0_mod = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
	for (int i = 0; i < 3; i++)iv[i] = rv0[i + 3] / v0_mod;
	for (int i = 0; i < 3; i++)dv[i] = dv_min * iv[i];

	tf = tmin;
	mf = m0 * exp(-fabs(dv_min * 1000.0) / g0 / Isp);

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
	//for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//升交和降交两支
		int bra = branch;
		NRi = NR;
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
			//if (as < Re_km + 200.0 || as > Re_km + 1000.0)continue;

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
			if (am < Re_km + 200.0) {
				flag0 = 0;
				return;
			}

			//此时am已求出，且i,Omega不变，需计算其他平轨道根数（近似解）
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			if (em < 0.0 || em > 1.0) {
				flag0 = 0;
				return;
			}
			if (am * (1.0 - em) < Re_km + 200.0) {
				flag0 = 0;
				return;
			}
			if (am * (1.0 + em) > Re_km + 1000.0) {
				flag0 = 0;
				return;
			}

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

			if ((double)sign * dv_mod < 0.0) {
				flag0 = 0;
				return;//新增脉冲方向的选择
			}

			//if (fabs(dv_min) > fabs(dv_mod))
			//if(fabs(tG_J2 - dt < 0.7 * 3600.0))
			dv_min = dv_mod;
			tmin = tG_J2;
			NR = NRi;
			flag0 = 1;
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
void single_imp_J2Lambert(double* dv, double* RVf, int& flag, const double* RV0, const double& t0, const double& tf, const double& lambda0, const double& phi0, const double& h, const int& N, const int& branch) {
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
	J2Lambert_short(flag_temp, RV1_temp, RVf_temp, RV0, tf - t0, mu_km_s, N, branch);

	bool ifvisible = is_target_visible(RVf_temp, R_Target, 19.5 * D2R);
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
//TODO: 修改优化时的propagate，使用J2线性模型
//优化变量（被摄动项）：
//		dv[3]：脉冲，每个分量扰动量30m/s
//		tf：实际的终末时刻，从6h开始扰动，区间为[0,12h]
//优化指标：
//		dv：脉冲大小
//参数：
//		t0：初始时刻
//		RV0[6]：初始位置速度
//		tf：优化前的终末时刻，固定3h
//		dv[3]：优化前的脉冲
//		target_id：要观测的目标序号（0-20）
//约束：
//		保证目标可见
//		轨道高度200-1000km
void perturbation(double* dv, double& tf, const std::vector<double>& X) {
	tf += (X[3] - 0.5) * 1.0e4;								
	dv[0] += (X[0] - 0.5) * 1.0e2;
	dv[1] += (X[1] - 0.5) * 1.0e2;
	dv[2] += (X[2] - 0.5) * 1.0e2;
}

//----------------------------------------------
// 帮助函数：给定 tf，计算对应的 angle
//----------------------------------------------
double calcAngle(int id,
	double tf,
	const double RV1[6],
	double t0)
{
	// 1. 计算目标点坐标
	double target_R[3];
	get_target_R(id, tf, target_R);

	// 2. 根据 RV1 传播到 tf 得到 RVf
	double RVf[6];
	propagate_linearJ2(RV1, RVf, t0, tf);

	// 3. 计算此时的视线角度
	//double angle = angle_target_visible(RVf, target_R, 0.0);
	double angle = angle_earth_view_angle(RVf, target_R);
	return angle;
}

//----------------------------------------------
// 黄金分割搜索核心函数：在 [ta, tb] 上搜索 angle 的最小值
//----------------------------------------------
double golden_section_search(int id,
	const double RV1[6],
	double t0,
	double ta,
	double tb,
	double tolerance = 1.0e-1)
{
	double angle_dv_min = 1e10;
	double min_angle = (ta + tb) / 2.0;
	//先粗筛
	double t_stp = 50.0;
	for (double x = ta; x <= tb; x += t_stp)
	{
		double angle = calcAngle(id, x, RV1, t0);
		if(angle < angle_dv_min )
		{
			angle_dv_min = angle;
			min_angle = x;
		}
	}

	ta = min_angle - t_stp;
	tb = min_angle + t_stp;

	// 黄金分割比
	const double phi = (1.0 + std::sqrt(5.0)) / 2.0;

	// 初始区间内的两个测试点
	double c = tb - (tb - ta) / phi;
	double d = ta + (tb - ta) / phi;

	// 当区间宽度大于设定的精度时，不断缩小区间
	while ((tb - ta) > tolerance)
	{
		double angleC = calcAngle(id, c, RV1, t0);
		double angleD = calcAngle(id, d, RV1, t0);

		// 如果在 c 处比 d 处的 angle 大，说明最优解在 [c, tb] 中
		// 否则在 [ta, d] 中
		if (angleC > angleD)
		{
			// [c, tb] 继续
			ta = c;
			c = d;
			d = ta + (tb - ta) / phi;
		}
		else
		{
			// [ta, d] 继续
			tb = d;
			d = c;
			c = tb - (tb - ta) / phi;
		}
	}

	// 返回区间中点作为近似最优解
	return 0.5 * (ta + tb);
}


//----------------------------------------------
// 【新增】函数：在 [ta, tb] 中做全局搜索
// 假设函数最多有两个局部极小值（两个谷），
// 则先粗采样，再对可疑区间用黄金分割找局部谷值，最后选全局最小
//----------------------------------------------
double global_min_search_two_valleys(
	int id,
	const double RV1[6],
	double t0,
	double ta,
	double tb,
	double coarseStep = 10.0,   // 粗采样步长，可根据需求调整
	double localTol = 1.0e-1  // 黄金分割搜索的容差
)
{
	// 1) 粗采样
	std::vector<double> xs;
	for (double x = ta; x <= tb; x += coarseStep)
	{
		xs.push_back(x);
	}
	if (xs.back() < tb)
	{
		xs.push_back(tb); // 保证把 tb 包含进去
	}

	// 2) 计算每个采样点的函数值
	std::vector<double> fs(xs.size());
	for (size_t i = 0; i < xs.size(); ++i)
	{
		fs[i] = calcAngle(id, xs[i], RV1, t0);
	}

	// 3) 根据采样结果，找可能存在“局部谷”的区间
	//    如果 f[i-1] >= f[i] <= f[i+1]，则在 [xs[i-1], xs[i+1]] 里可能有极小值
	//    注意边界情况可以单独处理
	std::vector<std::pair<double, double>> candidateIntervals;
	for (size_t i = 1; i + 1 < xs.size(); ++i)
	{
		double fPrev = fs[i - 1];
		double fCurr = fs[i];
		double fNext = fs[i + 1];

		// 判断 "先降后升" 或者 ">= <= >="
		if (( fCurr <= fPrev && fCurr <= fNext) || i ==1 || i+1 == xs.size()-1)
		{
			double left = xs[i - 1];
			double right = xs[i + 1];
			candidateIntervals.push_back({ left, right });
		}
	}

	// 如果采样步长比较大，或者函数的谷碰巧在边缘，也可能没被识别
	// 可酌情把整段 [ta, tb] 也加进来兜底
	if (candidateIntervals.empty())
	{
		candidateIntervals.push_back({ ta, tb });
	}

	// 4) 用黄金分割搜索在这些候选区间内找局部最小
	double bestX = ta;
	double bestF = std::numeric_limits<double>::infinity();

	for (auto& seg : candidateIntervals)
	{
		double segA = seg.first;
		double segB = seg.second;

		// 注意：若 segA==segB 或 segA>segB(极端情况)要做保护
		if (segA > segB) std::swap(segA, segB);
		if (std::fabs(segB - segA) < 1e-12)
		{
			// 区间太小就跳过或手动计算
			continue;
		}

		// 在这个区间做“单峰搜索”
		double xLocMin = golden_section_search(id, RV1, t0, segA, segB, localTol);
		double fLocMin = calcAngle(id, xLocMin, RV1, t0);

		// 如果比当前最优还小，则更新
		if (fLocMin < bestF)
		{
			bestF = fLocMin;
			bestX = xLocMin;
		}
	}

	// 返回 (x, f)
	return bestX;
}


double obj_func_shooting(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	double* para = static_cast<double*>(f_data);

	// 初值
	const double t0 = para[0];
	double tf = para[1];
	double dv[3] = { para[2], para[3], para[4] };
	const double RV0[6] = { para[5], para[6], para[7], para[8], para[9], para[10] };
	const int id = static_cast<int>(para[11]);

	perturbation(dv, tf, X);

	double penalty_time = 0.0;
	if (tf < t0) {
		//std::cout << "tf < t0，不符合要求" << std::endl;
		//std::cout << std::endl;
		penalty_time = t0 - tf;
		
	}

	double RVf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double RV1[6]; memcpy(RV1, RV0, 6 * sizeof(double));
	for (int i = 0; i < 3; i++) RV1[i + 3] += dv[i];

	double coe[6]; int flag;
	rv2coe(flag, coe, RV1, mu_km_s);
	double a = coe[0], e = coe[1];

	double impulse = V_Norm2(dv, 3);

	//ZZ Test
	//double geotic[2];
	//get_target_geogetic(id, tf, geotic);
	//modify_tg(t0, RV1, geotic[0], tf);


	// 用二体近似，为了给J2留出余量，设置20km的冗余
	double peri = a * (1 - e) - Re_km;
	double apo = a * (1 + e) - Re_km;
	double penalty_orbit = 0.0;
	penalty_orbit += fabs(std::min(peri, 200.0) - 200.0);
	penalty_orbit += fabs(std::max(apo, 1000.0) - 1000.0);


	//propagate_j2(RV1, RVf, t0, tf, 1e-10, 1e-10);
	//propagate_linearJ2(RV1, RVf, t0, tf);
	//rv02rvf(flag, RVf, RV1, tf - t0, mu_km_s);


	//ZZ Test
	// 调用黄金分割搜索，得到使 angle 最小的时间 tf_opt
	//tf = golden_section_search(id, RV1, t0, std::max (tf-1000.0, t0 + 60.0), tf + 1000.0, 1e-1);

	//tf = global_min_search_two_valleys(id, RV1, t0, std::max(tf - 1800.0, t0 + 60.0), tf + 1800.0, 100.0);

	//tf = 8200.0;

	propagate_linearJ2(RV1, RVf, t0, tf);

	double target_R[3];
	get_target_R(id, tf, target_R);
	//bool ifVisible = is_target_visible(RVf, target_R, 19.8 * D2R);
	double penalty_angle = 0.0;

	//if (!ifVisible) 
	{
		//double Rf[3] = { RVf[0], RVf[1], RVf[2] };
		//double Rf_norm = V_Norm2(Rf, 3);
		//double Rf_fixed[3];
		//V_Multi(Rf_fixed, Rf, Re_km / Rf_norm, 3);
		//double dR[3];
		//V_Minus(dR, Rf_fixed, target_R, 3);
		//impulse += V_Norm2(dR, 3) / Re_km / D2R;

		//penalty_angle = angle_target_visible(RVf, target_R, 19.8 * D2R);
		double Rf_norm = V_Norm2(RVf, 3);
		double cone_angle = 19.5 * D2R; //19.75
		double gamma_max = asin(sin(cone_angle) * Rf_norm / Re_km) - cone_angle;
		penalty_angle = std::max(angle_earth_view_angle(RVf, target_R) - gamma_max,0.0);
		impulse += penalty_angle * penalty_angle * 100000.0;
	}

	//para[1] = tf;

	return impulse + penalty_orbit* penalty_orbit * 100.0 +  penalty_time * penalty_time * 1.0e3;
}

//// TODO：优化之后，用高精度propagate再积一次分
//void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const int& target_id, const int& NR, const int& branch) {
//	flag = 0;
//	double target_geo[2];
//	double lambda, phi;  //lambda 经度，phi 纬度
//	if (target_id != 20) {
//		get_target_geogetic(target_id, tf, target_geo);
//		lambda = target_geo[1];
//		phi = target_geo[0];
//	}
//	else { // 船
//		get_target_geogetic(target_id, 2.0 * 86400.0, target_geo);
//		double phi_end = target_geo[0];
//		get_target_geogetic(target_id, 0.0, target_geo); //TODO： ZZ 计时是否错误？ 圈数应该是从开头到现在的圈数，或者用时间来计更合理吧？
//		double phi_start = target_geo[0];
//
//		// 纬度更精确一些：以一圈105min计，将纬度等分
//		double dphi = (phi_start - phi_end) * double(NR) * 6300.0 / 172800.0;
//		phi = phi_start + dphi;
//		lambda = target_geo[1];
//	}
//	
//	double m0 = 1000.0, mf; 
//
//	//注意张刚论文的方法里，tf是飞行时长，不是时刻
//	//对于地面目标，放宽打靶经纬度，先计算地心夹角范围，再根据这个给经度和纬度都撒网格，从中选脉冲最低的
//	//对于船，纬度已经比较近似
//	if (target_id != 20) {
//		// 计算地心角范围beta
//		double R0[3] = { RV0[0], RV0[1], RV0[2] };
//		double r = V_Norm2(R0, 3);
//		double gamma = 20.0 * D2R;
//		double beta = asin(r * sin(gamma) / Re_km) - gamma;
//		int mesh_size = 20, flag_temp = 0;
//		double dv_norm = 1.0e10;
//		double dv_temp[3], tf_temp = 0.0, lambda_temp, phi_temp;
//
//		single_imp(m0, t0, RV0, lambda, phi, 1, flag_temp, mf, tf_temp, dv_temp, NR, branch);
//		if (flag_temp == 1) {
//			if (V_Norm2(dv_temp, 3) < dv_norm) {
//				flag = 1;
//				dv_norm = V_Norm2(dv_temp, 3);
//				memcpy(dv, dv_temp, 3 * sizeof(double));
//				tf = tf_temp;
//			}
//		}
//
//		for (int i = 0; i < mesh_size; i++) {
//			lambda_temp = lambda + 2.0 * beta / mesh_size * (i - mesh_size / 2.0);
//
//			for (int j = 0; j < mesh_size; j++) {
//				phi_temp = phi + 2.0 * beta / mesh_size * (j - mesh_size / 2.0);
//				single_imp(m0, t0, RV0, lambda_temp, phi_temp, 1, flag_temp, mf, tf_temp, dv_temp, NR, branch);
//
//				if (flag_temp == 1) { 
//					if (V_Norm2(dv_temp, 3) < dv_norm) {
//						flag = 1;
//						dv_norm = V_Norm2(dv_temp, 3);
//						memcpy(dv, dv_temp, 3 * sizeof(double));
//						tf = tf_temp;
//					}
//				}
//			}
//		}
//	}
//	else { // 船
//		// 计算地心角范围beta
//		double R0[3] = { RV0[0], RV0[1], RV0[2] };
//		double r = V_Norm2(R0, 3);
//		double gamma = 20.0 * D2R;
//		double beta = asin(r * sin(gamma) / Re_km) - gamma;
//		int mesh_size = 20, flag_temp = 0;
//		double dv_norm = 1.0e10;
//		double dv_temp[3], tf_temp = 0.0, lambda_temp, phi_temp;
//		for (int i = 0; i < mesh_size; i++) {
//			lambda_temp = lambda + 2.0 * beta / mesh_size * (i - mesh_size / 2.0);
//
//			for (int j = 0; j < mesh_size; j++) {
//				phi_temp = phi + 2.0 * beta / mesh_size * (j - mesh_size / 2.0);
//				single_imp_ship(m0, t0, RV0, lambda_temp, phi_temp, 1, flag_temp, mf, tf_temp, dv_temp, NR, branch);
//
//				if (flag_temp == 1) {
//					if (V_Norm2(dv_temp, 3) < dv_norm) {
//						flag = 1;
//						dv_norm = V_Norm2(dv_temp, 3);
//						memcpy(dv, dv_temp, 3 * sizeof(double));
//						tf = tf_temp;
//					}
//				}
//			}
//		}
//	}
//	tf += t0;
//	
//	if (flag != 1) {
//		for (int j = 0; j < 3; j++) {
//			dv[j] = penalty;
//			RVf[j] = penalty; RVf[j + 3] = penalty;
//		}
//		tf = penalty;
//		return;
//	}
//
//	double f_data[12] = { t0, tf, dv[0], dv[1], dv[2], RV0[0], RV0[1], RV0[2], RV0[3], RV0[4], RV0[5], target_id};
//
//	double impulse = 0.0;
//	std::vector<double> X = { 0.5, 0.5, 0.5, 0.5 };
//	nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 10000);		//不输出
//
//	perturbation(dv, tf, X);
//
//	double RV1[6];
//	memcpy(RV1, RV0, 6 * sizeof(double));
//	for (int i = 0; i < 3; i++) RV1[i + 3] += dv[i];
//	propagate_j2(RV1, RVf, t0, tf);
//
//	double target_R[3];
//	get_target_R(target_id, tf, target_R);
//	bool ifVisible = is_target_visible(RVf, target_R, 20.0 * D2R);
//
//	if (!ifVisible) {
//		flag = 0;
//		for (int j = 0; j < 3; j++) {
//			dv[j] = penalty;
//			RVf[j] = penalty; RVf[j + 3] = penalty;
//		}
//		tf = penalty;
//		return;
//	}
//
//}


// TODO：优化之后，用高精度propagate再积一次分
void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const int& target_id, const int& NR, const int& branch) {
	flag = 0;
	double target_geo[2];
	double lambda, phi;  //lambda 经度，phi 纬度
	if (target_id != 20) {
		get_target_geogetic(target_id, tf, target_geo);
		lambda = target_geo[1];
		phi = target_geo[0];
	}
	else { // 船
		get_target_geogetic(target_id, 2.0 * 86400.0, target_geo);
		double phi_end = target_geo[0];
		get_target_geogetic(target_id, 0.0, target_geo); //TODO： ZZ 计时是否错误？ 圈数应该是从开头到现在的圈数，或者用时间来计更合理吧？
		double phi_start = target_geo[0];

		// 纬度更精确一些：以一圈105min计，将纬度等分
		double dphi = (phi_start - phi_end) * double(NR) * 6300.0 / 172800.0;
		phi = phi_start + dphi;
		lambda = target_geo[1];
	}

	double m0 = 1000.0, mf;

	//注意张刚论文的方法里，tf是飞行时长，不是时刻
	//对于地面目标，放宽打靶经纬度，先计算地心夹角范围，再根据这个给经度和纬度都撒网格，从中选脉冲最低的
	//对于船，纬度已经比较近似
	if (target_id != 20) {
		// 计算地心角范围beta
		double R0[3] = { RV0[0], RV0[1], RV0[2] };
		double r = V_Norm2(R0, 3);
		double gamma = 35.0 * D2R;
		double beta = asin(r * sin(gamma) / Re_km) - gamma;
		int mesh_size = 20, flag_temp = 0;
		double dv_norm = 1.0e10;
		double dv_temp[3], tf_temp = 0.0, lambda_temp, phi_temp;

		//ZZ test
		//single_imp(m0, t0, RV0, lambda, phi, 1, flag_temp, mf, tf_temp, dv_temp, NR, branch);
		single_imp_zzmodified(m0, t0, RV0, lambda, phi, 1, flag_temp, mf, tf_temp, dv_temp, branch);

		if (flag_temp == 1) {
			if (V_Norm2(dv_temp, 3) < dv_norm) {
				flag = 1;
				dv_norm = V_Norm2(dv_temp, 3);
				memcpy(dv, dv_temp, 3 * sizeof(double));
				tf = tf_temp;
			}
		}

		//ZZ test

		//for (int i = 0; i < mesh_size; i++) {
		//	lambda_temp = lambda + 2.0 * beta / mesh_size * (i - mesh_size / 2.0);

		//	for (int j = 0; j < mesh_size; j++) {
		//		phi_temp = phi + 2.0 * beta / mesh_size * (j - mesh_size / 2.0);
		//		//single_imp(m0, t0, RV0, lambda_temp, phi_temp, 1, flag_temp, mf, tf_temp, dv_temp, NR, branch);
		//		single_imp_zzmodified(m0, t0, RV0, lambda_temp, phi_temp, 1, flag_temp, mf, tf_temp, dv_temp,  branch);
		//		if (flag_temp == 1) {
		//			if (V_Norm2(dv_temp, 3) < dv_norm) {
		//				flag = 1;
		//				dv_norm = V_Norm2(dv_temp, 3);
		//				memcpy(dv, dv_temp, 3 * sizeof(double));
		//				tf = tf_temp;
		//			}
		//		}
		//	}
		//}
	}
	else { // 船
		// 计算地心角范围beta
		double R0[3] = { RV0[0], RV0[1], RV0[2] };
		double r = V_Norm2(R0, 3);
		double gamma = 20.0 * D2R;
		double beta = asin(r * sin(gamma) / Re_km) - gamma;
		int mesh_size = 20, flag_temp = 0;
		double dv_norm = 1.0e10;
		double dv_temp[3], tf_temp = 0.0, lambda_temp, phi_temp;
		//ZZ test
		//for (int i = 0; i < mesh_size; i++) {
		//	lambda_temp = lambda + 2.0 * beta / mesh_size * (i - mesh_size / 2.0);

		//	for (int j = 0; j < mesh_size; j++) {
		//		phi_temp = phi + 2.0 * beta / mesh_size * (j - mesh_size / 2.0);
		//		single_imp_ship(m0, t0, RV0, lambda_temp, phi_temp, 1, flag_temp, mf, tf_temp, dv_temp, NR, branch);

		//		if (flag_temp == 1) {
		//			if (V_Norm2(dv_temp, 3) < dv_norm) {
		//				flag = 1;
		//				dv_norm = V_Norm2(dv_temp, 3);
		//				memcpy(dv, dv_temp, 3 * sizeof(double));
		//				tf = tf_temp;
		//			}
		//		}
		//	}
		//}
	}
	tf += t0;

	if (flag != 1) {
		for (int j = 0; j < 3; j++) {
			dv[j] = penalty;
			RVf[j] = penalty; RVf[j + 3] = penalty;
		}
		tf = penalty;
		return;
	}

	double f_data[12] = { t0, tf, dv[0], dv[1], dv[2], RV0[0], RV0[1], RV0[2], RV0[3], RV0[4], RV0[5], target_id };

	double impulse = 0.0;
	std::vector<double> X = { 0.5, 0.5, 0.5, 0.5};
	nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 100000);		//不输出

	//PSO(obj_func_shooting, f_data, X, impulse, X.size(), 0, 100000, 2000000);		//不输出



	perturbation(dv, tf, X);

	double RV0_temp[6]; memcpy(RV0_temp, RV0, 6 * sizeof(double));
	for (int i = 0; i < 3; i++) RV0_temp[i + 3] += dv[i];
	//modify_tg(t0, RV0_temp, phi, tf);

	//tf = golden_section_search(target_id, RV0_temp, t0, std::max(tf - 500.0, t0 + 60.0), tf + 500.0, 1e-1);


	double RV1[6];
	memcpy(RV1, RV0, 6 * sizeof(double));
	for (int i = 0; i < 3; i++) RV1[i + 3] += dv[i];
	propagate_j2(RV1, RVf, t0, tf);

	double target_R[3];
	get_target_R(target_id, tf, target_R);
	bool ifVisible = is_target_visible(RVf, target_R, 20.0 * D2R);

	if (!ifVisible) {
		flag = 0;
		for (int j = 0; j < 3; j++) {
			dv[j] = penalty;
			RVf[j] = penalty; RVf[j + 3] = penalty;
		}
		tf = penalty;
		return;
	}

}


std::vector<dv_shooding_info> obs_shooting_zzmodified(const double& t0, const double* RV0, const int& target_id) {

	int flag = 0;
	double target_geo[2];
	double lambda, phi;  //lambda 经度，phi 纬度
	if (target_id != 20) {
		get_target_geogetic(target_id, t0, target_geo);
		lambda = target_geo[1];
		phi = target_geo[0];
	}
	else { // 船
		get_target_geogetic(target_id, t0, target_geo);
		lambda = target_geo[1];
		phi = target_geo[0];
	}

	//double m0 = 1000.0, mf;
	//double tg = 0.0; // 注意是经过的时间，不是时刻
	//double dv[3];
	//注意tf是飞行时长，不是时刻
	// {
	//	// 计算地心角范围beta
	//	double R0[3] = { RV0[0], RV0[1], RV0[2] };
	//	double r = V_Norm2(R0, 3);
	//	int  flag_temp = 0;
	//	double dv_norm = 1.0e10;
	//	double dv_temp[3];

	//	//ZZ test
	//	single_imp_zzmodified(m0, t0, RV0, lambda, phi, 1, flag_temp, mf, tg, dv_temp, branch);

	//	if (flag_temp == 1) {
	//		if (V_Norm2(dv_temp, 3) < dv_norm) {
	//			flag = 1;
	//			dv_norm = V_Norm2(dv_temp, 3);
	//			memcpy(dv, dv_temp, 3 * sizeof(double));
	//		}
	//	}
	//}

	//double tf = t0 + tg;

	double tf = t0 + 5400;

	std::vector<dv_shooding_info> dv_infos;
	//if (flag != 1) {return dv_infos;}

	double step_t = 900.0;

	for (int i = 0; i < 11; i++)
	{
		dv_shooding_info dv_info_temp;
		dv_info_temp.tf = tf + (i-5)* step_t;
		int flag_lambert;
		lambert_dv(flag_lambert, t0, RV0,  dv_info_temp.dv, dv_info_temp.tf,  target_id);

		//memcpy(dv_info_temp.dv, dv, 3 * sizeof(double));
		if(dv_info_temp.tf < t0 + 30.0 ) continue;

		double f_data[12] = { t0, dv_info_temp.tf, dv_info_temp.dv[0], dv_info_temp.dv[1], dv_info_temp.dv[2], RV0[0], RV0[1], RV0[2], RV0[3], RV0[4], RV0[5], (double)target_id };
		double impulse = 0.0;
		std::vector<double> X = { 0.5, 0.5, 0.5, 0.5 };
		nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 1000);		//不输出
		perturbation(dv_info_temp.dv, dv_info_temp.tf, X);

		//验证是否能可见，且注意要整分钟下
		double RV0_temp[6];
		memcpy(RV0_temp, RV0, 6 * sizeof(double));
		for (int i = 0; i < 3; i++) RV0_temp[i + 3] += dv_info_temp.dv[i];
		double tf_60 = round(dv_info_temp.tf / 60) * 60.0;  //取整到分钟，跟问题设定一致
		dv_info_temp.tf = tf_60;
		propagate_j2(RV0_temp, dv_info_temp.RVF, t0, dv_info_temp.tf);

		double target_R[3];
		get_target_R(target_id, dv_info_temp.tf, target_R);
		bool ifVisible = is_target_visible(dv_info_temp.RVF, target_R, 19.99 * D2R);

		if (!ifVisible) { continue; }

		dv_infos.push_back(dv_info_temp);
	}
	return  dv_infos;
}