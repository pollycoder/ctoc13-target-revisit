#include "single_impluse.h"
extern int solve_cubic(double a, double b, double c, double d, double* xout);


//���������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		lambda:		Ŀ���ĵ��澭�ȣ�rad
//		phi:		Ŀ���ĵ���γ�ȣ�rad
//		Day:		�趨�����������㷶Χ(Day-1)~Day
//		branch:		������0�򽵽���1
//		sign:		��Ҫ��dv�ķ��ţ�1����-1��
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
				int& flag0, double& mf, double& tf, double* dv, int& NR,const int branch, const int sign)
{
	double lambda, phi;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6],ocoe0[6];//��ʼ�������
	double c1, c2, k1, k2, k3, A, B, C, D;//�м����
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//������ز�����J2�¶�ʱ���һ�׵���
	double mcoe[6], ocoe[6];//���ս����ƽ˲��ת��
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//����J2����ʵ��λ���ٶ�

	flag0 = 0;

	//rv2coe(flag, ocoe0, rv0, mu_km_s);//�������ת���ɹ������
	//ocoe0[0] *= 1000.0;
	//O2M(ocoe0, coe0, J2);
	//coe0[0] /= 1000.0;
	rv2coe(flag, coe0, rv0, mu_km_s);//�������ת���ɹ������

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);
	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	//lambda = lambda0 * D2R;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	//phi = phi0 * D2R;
	lambda = lambda0;
	phi = phi0;

	Nmax = floor(double(Day) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//���ת��Ȧ��
	Nmin = floor((double(Day) - 1.0) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//����ת��Ȧ��

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//�����ͽ�����֧
		int bra = branch;
		{
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//������(����as)�����η��̣�����fabs(as-a0)��С������
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

			//����J2ģ�����
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);

			//����J2ģ���µ����η������am����ƽ�볤��
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

			//��ʱam���������i,Omega���䣬���������ƽ������������ƽ⣩
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			//if (em < 0.0 || em > 1.0)continue;
			if (am * (1.0 - em) < Re_km + 200.0)continue;
			//if (am * (1.0 + em) > Re_km + 1000.0)continue;

			//ƽ��ת��Ϊ˲��������ƽ˲��ת��������ʹ�õ�ReΪm��������Ҳ��ʱ�ĳ�m������
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			//�����ʼʱ�̵���������
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0]/1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

			if ((double)sign * dv_mod < 0.0)continue;//�������巽���ѡ��

			if (fabs(dv_min) > fabs(dv_mod))
			{
				dv_min = dv_mod;

				//C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(am, -3.5);
				//Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - em * em) / (1.0 - em * em);
				//tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);//��ȷ��J2����ʱ��
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


//���������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		lambda:		Ŀ���ĵ��澭�ȣ�rad
//		phi:		Ŀ���ĵ���γ�ȣ�rad
//		Day:		�趨�����������㷶Χ(Day-1)~Day
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR) 
{
	double lambda, phi;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6], ocoe0[6];//��ʼ�������
	double c1, c2, k1, k2, k3, A, B, C, D;//�м����
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//������ز�����J2�¶�ʱ���һ�׵���
	double mcoe[6], ocoe[6];//���ս����ƽ˲��ת��
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//����J2����ʵ��λ���ٶ�

	flag0 = 0;

	//rv2coe(flag, ocoe0, rv0, mu_km_s);//�������ת���ɹ������
	//ocoe0[0] *= 1000.0;
	//O2M(ocoe0, coe0, J2);
	//coe0[0] /= 1000.0;
	rv2coe(flag, coe0, rv0, mu_km_s);//�������ת���ɹ������

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);
	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	//lambda = lambda0 * D2R;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	//phi = phi0 * D2R;
	lambda = lambda0;
	phi = phi0;

	Nmax = floor(double(Day) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//���ת��Ȧ��
	Nmin = floor((double(Day) - 1.0) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//����ת��Ȧ��

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		for (int bra = 0; bra < 2; bra++)//�����ͽ�����֧
		//int bra = branch;
		{
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//������(����as)�����η��̣�����fabs(as-a0)��С������
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

			//����J2ģ�����
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);

			//����J2ģ���µ����η������am����ƽ�볤��
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

			//��ʱam���������i,Omega���䣬���������ƽ������������ƽ⣩
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			//if (em < 0.0 || em > 1.0)continue;
			if (am * (1.0 - em) < Re_km + 200.0)continue;
			if (am * (1.0 + em) > Re_km + 1500.0)continue;

			//ƽ��ת��Ϊ˲��������ƽ˲��ת��������ʹ�õ�ReΪm��������Ҳ��ʱ�ĳ�m������
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			//�����ʼʱ�̵���������
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

			//if ((double)sign * dv_mod < 0.0)continue;//�������巽���ѡ��

			if (fabs(dv_min) > fabs(dv_mod))
			{
				dv_min = dv_mod;

				//C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(am, -3.5);
				//Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - em * em) / (1.0 - em * em);
				//tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);//��ȷ��J2����ʱ��
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
	double lambda, phi;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6], ocoe0[6];//��ʼ�������
	double c1, c2, k1, k2, k3, A, B, C, D;//�м����
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//������ز�����J2�¶�ʱ���һ�׵���
	double mcoe[6], ocoe[6];//���ս����ƽ˲��ת��
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//����J2����ʵ��λ���ٶ�

	flag0 = 0;

	//rv2coe(flag, ocoe0, rv0, mu_km_s);//�������ת���ɹ������
	//ocoe0[0] *= 1000.0;
	//O2M(ocoe0, coe0, J2);
	//coe0[0] /= 1000.0;
	rv2coe(flag, coe0, rv0, mu_km_s);//�������ת���ɹ������

	a0 = coe0[0];
	e0 = coe0[1];
	inc0 = coe0[2];
	Omega0 = coe0[3];
	omega0 = coe0[4];
	f00 = coe0[5];
	E00 = f2E(flag, f00, e0);
	M00 = E2M(flag, E00, e0);
	gamma00 = atan(e0 * sin(f00) / (1.0 + e0 * cos(f00)));

	//lambda = lambda0 * D2R;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	//phi = phi0 * D2R;
	lambda = lambda0;
	phi = phi0;

	Nmax = floor(dt / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s))) + 1;//���ת��Ȧ��
	Nmin = floor(dt / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//����ת��Ȧ��

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//�����ͽ�����֧
		int bra = branch;
		{
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_)) / omegaE;

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//������(����as)�����η��̣�����fabs(as-a0)��С������
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

			//����J2ģ�����
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_) / (omegaE - Omega_J2_);

			//����J2ģ���µ����η������am����ƽ�볤��
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

			//��ʱam���������i,Omega���䣬���������ƽ������������ƽ⣩
			em = e0 + k2 * (am / a0 - 1.0);
			omegam = omega0 + k1 * (am / a0 - 1.0);
			M0m = M00 - k3 * (am / a0 - 1.0);
			//if (em < 0.0 || em > 1.0)continue;
			//if (am * (1.0 - em) < Re_km + 200.0)continue;
			//if (am * (1.0 + em) > Re_km + 1000.0)continue;

			//ƽ��ת��Ϊ˲��������ƽ˲��ת��������ʹ�õ�ReΪm��������Ҳ��ʱ�ĳ�m������
			mcoe[0] = am * 1000.0;
			mcoe[1] = em;
			mcoe[2] = inc0;
			mcoe[3] = Omega0;
			mcoe[4] = omegam;
			mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
			M2O(mcoe, ocoe, J2);

			//�����ʼʱ�̵���������
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

			if ((double)sign * dv_mod < 0.0)continue;//�������巽���ѡ��

			if (fabs(dv_min) > fabs(dv_mod))
			//if(fabs(tG_J2 - dt < 0.7 * 3600.0))
			{
				dv_min = dv_mod;

				//C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(am, -3.5);
				//Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - em * em) / (1.0 - em * em);
				//tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);//��ȷ��J2����ʱ��
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
	tf += (X[3] - 0.5) * 2 * hour;															// Time: ��0.5h

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


//CTOC13������J2Lambert���⣬���е�����������
//���Ÿ������������ƣ��������岢���ǽ���������һ������߶ȷ�Χ�ڽ���ѡ������������ܹ۲⵽����Ŀ���
//���룺
//		RV0[6]����ʼλ���ٶ�
//		t0����ʼʱ��
//		tf���ն�ʱ��
//�����
//		dv[3]�����յ�����
//		RVf[6]���ն�λ���ٶ�
//		flag�����ɹ�����1�����ʧ�ܷ���0    
//		h�����ղ��õĸ߶�   
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

	//������ɹ���ͷ�
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

	//����ɹ����������
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

//CTOC13������J2Lambert���Ż�����Ŀ���Ĺ���߶�
//�Ż����������㶯���
//		h������߶�
//		(lambda, phi)��ʵ�ʾ����ĵ��澭γ�ȣ����Ŷ�ʱ�䣬�ٽ�һ���Ŷ���γ�ȣ��Ŷ���180��
//		tf��ʵ�ʵ���ĩʱ�̣���3h��ʼ�Ŷ�������Ϊ[0,6h]
//�Ż�ָ�꣺
//		dv�������С
//������
//		t0����ʼʱ��
//		RV0[6]����ʼλ���ٶ�
//		tf���Ż�ǰ����ĩʱ��
//		h0���Ż�ǰ�Ĺ���߶�
//		target_id��Ҫ�۲��Ŀ����ţ�0-20��
//Լ����
//		��֤Ŀ��ɼ�
//		����߶�200-1000km
//��ʽ��PSO��Ӧ��Obj_func��ȫһ�£�����ֻ��Ҫ��ȡ���յ��Ż�����X
//��ͬ��֮ǰÿ�ζ�Ҫ��дget_value�ٰ�װ������Ϊ��childnode����������˼·������
void perturbation(double& h, double& tf, double& lambda, double& phi, const std::vector<double>& X, const double& h0, const int& id) {
	const double double_h_pert = 100.0;
	const double double_angle_pert = 12.0 * D2R;
	h = h0 + (X[0] - 0.5) * double_h_pert;
	tf += (X[3] - 0.5) * 43200.0;								// �Ŷ���γ��֮ǰ�����Ŷ�ʱ��

	// ��ȡʱ���Ŷ���ľ�γ�ȣ��Ե���Ŀ����Ӱ�죬��ҪӰ�캣��Ŀ�꣩
	double target_geo[2];
	get_target_geogetic(id, tf, target_geo);
	const double lambda0 = target_geo[1];
	const double phi0 = target_geo[0];

	// ��ʱ���Ŷ���ľ�γ�������Ƕ��Ŷ���������Ŀ�궼��Ӱ�죩
	lambda = lambda0 + (X[1] - 0.5) * double_angle_pert;
	phi = phi0 + (X[2] - 0.5) * double_angle_pert;
}


double obj_func_shooting(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	double* para = static_cast<double*>(f_data);

	// ��ֵ
	const double t0 = para[0];
	double tf = para[1];
	const double h0 = para[2];
	const double RV0[6] = { para[3], para[4], para[5], para[6], para[7], para[8] };
	const int id = static_cast<int>(para[9]);

	double lambda, phi, h;
	perturbation(h, tf, lambda, phi, X, h0, id);

	if (tf < t0) {
		//std::cout << "tf < t0��������Ҫ��" << std::endl;
		//std::cout << std::endl;
		return penalty;
	}

	double dv[3] = { 0.0, 0.0, 0.0 };
	double impulse = penalty;
	double RVf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	int flag = -1;
	single_imp(dv, RVf, flag, RV0, t0, tf, lambda, phi, h);

	if(h < 220.0 || h > 980.0) { 
		//std::cout << "�߶Ȳ�����Ҫ��, h = " << h << " km" << std::endl;
		//std::cout << std::endl;
		return penalty; 
	}
	if(flag != 1) { 
		//std::cout << "J2 Lambert ���ʧ��" << std::endl;
		//std::cout << std::endl;
		return penalty; 
	}

	double target_R[3];
	get_target_R(id, tf, target_R);
	bool ifVisible = is_target_visible(RVf, target_R, 19.5 * D2R);
	if(!ifVisible) { 
		//std::cout << "��к�Ŀ�겻�ɼ�" << std::endl;
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
	nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 100);		//�����

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