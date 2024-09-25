#include<math.h>
#include"single_impluse.h"
#include"Constant.h"
#include"OrbitFun.h"
#include"propagate.h"
#include"OrbitMath.h"

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
			if (as < Re_km || as > Re_km + 1000.0)continue;

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
			if (am * (1.0 + em) > Re_km + 1000.0)continue;

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


//��Ŀ��ĵ�����ת���ٶȼӳɣ�rad/s
const double omegaE_plus[5] = { (m_long[0][1] - m_long[0][0]) / oneway[0] * D2R, (m_long[1][1] - m_long[1][0]) / oneway[1] * D2R,
								(m_long[2][1] - m_long[2][0]) / oneway[2] * D2R, (m_long[3][1] - m_long[3][0]) / oneway[3] * D2R,
								(m_long[4][1] - m_long[4][0]) / oneway[4] * D2R };


//����ʵʱ��Ŀ��λ�ã�rad
int GetMov(double t, int ID, double& lat, double& lon)
{
	const double boundary[5][3] = { 37.7000, 140.9629, -122.6165 + 360, -38.7300, -73.3911 + 360, 177.9531, 40.6000,  -73.2332 + 360, -9.2553 + 360,
						 -22.6866, 14.3732 + 360, -42.0198 + 360, -6.3085,  105.2459, 39.0247 };//γ�ȡ����ȡ����ȣ������Ȼ��㵽[30, 390]
	const double t_ow[5] = { 8.0 * 86400.0, 8.0 * 86400.0, 5.0 * 86400.0, 6.0 * 86400.0, 7.0 * 86400.0 };//����ʱ��

	if (t < 180.0 * 86400.0 || t > 240.0 * 86400.0)return 0;

	if (ID >= 201 && ID <= 205)
	{
		int id = ID - 201;
		double t_, start, end;

		t_ = t - 180.0 * 86400.0;
		while (t_ > 2 * t_ow[id])t_ -= 2 * t_ow[id];//ʱ�任��
		if (t_ < t_ow[id])//�ж��˶��������յ�
		{
			start = boundary[id][1];
			end = boundary[id][2];
		}
		else
		{
			t_ -= t_ow[id]; 
			start = boundary[id][2];
			end = boundary[id][1];
		}

		lat = boundary[id][0] * D2R;
		lon = start + (end - start) * t_ / t_ow[id];
		if (lon > 180.0)lon -= 360.0;
		lon *= D2R;

		return 1;
	}
	else
		return 0;
}


//��tʱ�̵�ID��Ŀ�갴�յ�Ч��ת���ٶ�rad/s��Ӧ��0ʱ�̣��������Ч��ֹ����ʼʱ�̣���γ��rad
int GetStart(int ID, double t, double& omegaE_eq, double& alpha_G0_eq, double& lambda, double& phi)
{
	if (t < 180.0 * 86400.0)return 0;
	if (ID < 201 || ID > 205)return 0;

	int mID, way_num;

	mID = ID - 201;
	way_num = floor((t - 180.0 * 86400.0) / oneway[mID]);

	phi = m_lat[mID] * D2R;
	if (way_num % 2 == 0)//����A���ӳ��
	{
		omegaE_eq = omegaE + omegaE_plus[mID];
		alpha_G0_eq = alpha_G0 + omegaE * (oneway[mID] * double(way_num) + 180.0 * 86400.0) - omegaE_eq * (oneway[mID] * double(way_num) + 180.0 * 86400.0);//��ʼGMST�ĵ�Ч
		while (alpha_G0_eq > D2PI)alpha_G0_eq -= D2PI;
		while (alpha_G0_eq < 0)alpha_G0_eq += D2PI;

		lambda = m_long[mID][0] * D2R;// -omegaE_eq * oneway[mID] * double(way_num);//��ʼ���ȵĵ�Ч
		while (lambda < -DPI)lambda += D2PI;
	}
	else//����
	{
		omegaE_eq = omegaE - omegaE_plus[mID];
		alpha_G0_eq = alpha_G0 + omegaE * (oneway[mID] * double(way_num + 1) + 180.0 * 86400.0) - omegaE_eq * (oneway[mID] * double(way_num + 1) + 180.0 * 86400.0);//��ʼGMST�ĵ�Ч
		while (alpha_G0_eq > D2PI)alpha_G0_eq -= D2PI;
		while (alpha_G0_eq < 0)alpha_G0_eq += D2PI;

		lambda = m_long[mID][0] * D2R;
		while (lambda < -DPI)lambda += D2PI;
	}


	return 1;
}

//���ٶȼӳɵķ���
//�ƶ�Ŀ��ĵ��������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		ID:			��Ŀ���ţ�201~205
//		Day:		�趨��������Ϊ��֤����ͨ��ȡ1���2��
//		branch:		������0�򽵽���1
//�����
//		flag:		����ɹ�����1��ʧ�ܷ���0����������ʱ�䷵��-1
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
//		rvt[6]		�������������NULL����
void single_imp2(const double m0, const double t0, const double* rv0, const int ID, const int Day,
				int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, double* rvt)
{
	int mID, way_num;
	double lambda, phi, omegaE_eq, alpha_G0_eq, m_t0;

	//ʱ�䲻������180��
	flag0 = 0;
	if (t0 < 180.0 * 86400.0)return;
	////�ƶ�Ŀ��λ��ȡ�õ��̵���㾭γ�ȣ��ƶ��ٶȵ�Ч��������ת��
	mID = ID - 201;
	way_num = floor((t0 - 180.0 * 86400.0) / oneway[mID]);
	//m_t0 = double(way_num) * oneway[mID];
	//phi = m_lat[mID] * D2R;
	//if (way_num % 2 == 0)//����
	//{
	//	lambda = m_long[mID][0] * D2R;
	//	omegaE_eq = omegaE + omegaE_plus[mID];
	//}
	//else//����
	//{
	//	lambda = m_long[mID][1] * D2R;
	//	omegaE_eq = omegaE - omegaE_plus[mID];
	//}
	if (!GetStart(ID, t0, omegaE_eq, alpha_G0_eq, lambda, phi))return;


	//������������
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6];//��ʼ�������
	double c1, c2, k1, k2, k3, A, B, C, D;//�м����
	double tG, tG_J2, alphaG0_, alphaGt, alphaGf;//GMST
	double C_J2, fs, as, am, delta_a, omega_ft, em, omegam, M0m, p0, dv_mod, v0_mod, iv[3];
	int flag, root_num, NRi, Nmax, Nmin;
	double dv_min, tmin;
	double root[3];
	double omega_J2_, Omega_J2_, M_J2_;//������ز�����J2�¶�ʱ���һ�׵���
	double mcoe[6], ocoe[6];//���ս����ƽ˲��ת��
	double rv0_[6], rvf[6], hf, rtar_F[3], rtar[3];//����J2����ʵ��λ���ٶ�

	int bra = branch;

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

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));


	Nmax = floor(double(Day) * D2PI / omegaE_eq / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//���ת��Ȧ��
	Nmin = floor((double(Day) - 1.0) * D2PI / omegaE_eq / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//����ת��Ȧ��

	alphaG0_ = alpha_G0_eq + omegaE_eq * t0;//��ʼʱ��GMST
	//alphaG0_ = alpha_G0 + omegaE * (180.0 * 86400.0 + m_t0) + omegaE_eq * (t0 - 180.0 * 86400.0 - m_t0);
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;
	omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

	c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
	c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
	alphaGt = Omega0 + atan2(c1, c2);
	//tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;
	while (alphaGt < alphaG0_)alphaGt += D2PI;
	while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
	tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE_eq;

	if (tG > double(way_num + 1.0) * oneway[mID])//ĩ�˳����������ڵ����
	{
		if (!GetStart(ID, t0 + tG, omegaE_eq, alpha_G0_eq, lambda, phi))return;//����ĩ��ʱ���������µĳ�ʼ��γ��

		Nmax = floor(double(Day) * D2PI / omegaE_eq / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//���ת��Ȧ��
		Nmin = floor((double(Day) - 1.0) * D2PI / omegaE_eq / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//����ת��Ȧ��

		alphaG0_ = alpha_G0_eq + omegaE_eq * t0;//��ʼʱ��GMST
		while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

		fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;
		omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

		c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
		c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
		alphaGt = Omega0 + atan2(c1, c2);
		while (alphaGt < alphaG0_)alphaGt += D2PI;
		while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
		tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE_eq;

		if (tG < double(way_num + 1.0) * oneway[mID])
		{
			flag0 = -1;
			return;//�Է�����
		}
	}


	dv_min = 1.0e10;
	tmin = 0.0;
	//��Ȧ����
	for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
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
		if (as < Re_km)continue;

		//����J2ģ�����
		C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
		omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
		Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
		M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
		//tG_J2 = (Omega0 + atan2(c1, c2) - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE - Omega_J2_);
		tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE_eq - Omega_J2_);


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
		if (fabs(em) > 0.1)continue;
		if (am * (1.0 - em) < Re_km + 200.0)continue;
		if (am * (1.0 + em) > Re_km + 1000.0)continue;

		//ƽ��ת��Ϊ˲��������ƽ˲��ת��������ʹ�õ�ReΪm��������Ҳ��ʱ�ĳ�m������
		mcoe[0] = am * 1000.0;
		mcoe[1] = em;
		mcoe[2] = inc0;
		mcoe[3] = Omega0;
		mcoe[4] = omegam;
		mcoe[5] = E2f(flag, M2E(flag, M0m, em, 100, 1.0e-14), em);
		M2O(mcoe, ocoe, J2);
		if (ocoe[0] * (1.0 - ocoe[1]) / 1000.0 < Re_km + 220.0)continue;//200km���Ǻ�Σ�յ�

		//�����ʼʱ�̵���������
		p0 = a0 * (1 - e0 * e0);
		dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));

		//if (dv_mod < -0.02)continue;//����Ҳ�Ƚ�Σ��

		if (fabs(dv_min) > fabs(dv_mod) && fabs(dv_mod) < 0.12)
		{
			//����ĩ��ʱ�̵�����״̬���ж�ʵ�ʸ߶�
			v0_mod = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
			for (int i = 0; i < 3; i++)iv[i] = rv0[i + 3] / v0_mod;
			for (int i = 0; i < 3; i++)dv[i] = dv_mod * iv[i];
			for (int i = 0; i < 3; i++)
			{
				rv0_[i] = rv0[i] * 1000.0;
				rv0_[i + 3] = (rv0[i + 3] + dv[i]) * 1000.0;
			}
			j2rv02rvf(rv0_, tG_J2, rvf);
			for (int i = 0; i < 6; i++)rvf[i] /= 1000.0;
			hf = sqrt(rvf[0] * rvf[0] + rvf[1] * rvf[1] + rvf[2] * rvf[2]) - Re_km;//ʵ�ʹ���߶ȣ���λ��km

			if (hf > 500.0 || hf < 200.0)
				continue;
			else
			{
				dv_min = dv_mod;
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
	{
		return;
	}
	else
	{
		////����ĩ��ʱ�̵�����״̬��ʵ�ʡ�����
		//for (int i = 0; i < 3; i++)
		//{
		//	rv0_[i] = rv0[i] * 1000.0;
		//	rv0_[i + 3] = (rv0[i + 3] + dv[i]) * 1000.0;
		//}
		//j2rv02rvf(rv0_, tf, rvf);
		//for (int i = 0; i < 6; i++)rvf[i] /= 1000.0;
		//hf = sqrt(rvf[0] * rvf[0] + rvf[1] * rvf[1] + rvf[2] * rvf[2]) - Re_km;//ʵ�ʹ���߶ȣ���λ��km

		//if (hf > 500.0 || hf < 200.0)return;

		//rtar_F[0] = Re_km * cos(phi) * cos(lambda);//����Ŀ���λ������ת��
		//rtar_F[1] = Re_km * cos(phi) * sin(lambda);
		//rtar_F[2] = Re_km * sin(phi);
		////alphaGf = alphaG0_ + omegaE * tf;//tfʱ��GMST
		//alphaGf = alphaG0_ + omegaE_eq * tf;

		//rtar[0] = rtar_F[0] * cos(alphaGf) - rtar_F[1] * sin(alphaGf);
		//rtar[1] = rtar_F[0] * sin(alphaGf) + rtar_F[1] * cos(alphaGf);
		//rtar[2] = rtar_F[2];
		//for (int i = 0; i < 3; i++)
		//{
		//	rvt[i] = rtar[i] * (1.0 + hf / Re_km);
		//	rvt[i + 3] = rvf[i + 3];
		//}

		flag0 = 1;

	}

};


//��������ת��ʱ�䣬���ݶ��忪����ģ��
//	ID:			�ƶ�Ŀ��201~205
//	t0:			��ʼʱ��
//	coe0[6]:	��ʼ�������
//	day:		ת��������һ��ȡ1
//	branch:		��֧
double tG_estimate(int ID, double t0, const double* coe0, int day, int branch, int Maxiter=10, double tol=1e-6)
{
	double phi, lambda, omega_ft, Omega0, inc0, alphaGt, alphaG0_, tG;
	double x1, x2, y1, y2, xnew, ynew, err, c1, c2;
	int iter, bra, flag;

	Omega0 = coe0[3];
	inc0 = coe0[2];
	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	bra = branch;
	//for (bra = 0; bra < 2; bra++)
	{
		x1 = 0.1 * 86400.0;//����ֵ
		flag = GetMov(t0 + x1, ID, phi, lambda);
		omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��
		c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
		c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
		alphaGt = Omega0 + atan2(c1, c2);
		while (alphaGt < alphaG0_)alphaGt += D2PI;
		while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
		y1 = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(day) - 1.0)) - x1 * omegaE;

		x2 = double(day) * 86400.0;//��������
		flag = GetMov(t0 + x2, ID, phi, lambda);
		omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��
		c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
		c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
		alphaGt = Omega0 + atan2(c1, c2);
		while (alphaGt < alphaG0_)alphaGt += D2PI;
		while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
		y2 = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(day) - 1.0)) - x2 * omegaE;

		//Newton iterations
		err = 1.0;
		iter = 0;
		while ((err > tol) && (y1 != y2) && (iter <= Maxiter))
		{
			iter++;
			xnew = (x1 * y2 - y1 * x2) / (y2 - y1);
			flag = GetMov(t0 + xnew, ID, phi, lambda);
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��
			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			ynew = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(day) - 1.0)) - xnew * omegaE;

			x1 = x2;
			y1 = y2;
			x2 = xnew;
			y2 = ynew;
			err = fabs(x1 - xnew);
		}
		tG = xnew;

		//if (tG > 0.5 * 86400.0)return tG;
	}

	return tG;
}

//��δ��ɵ��ԣ���������tG���ƶ�Ŀ��λ��
//�ƶ�Ŀ��ĵ��������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		ID:			�ƶ�Ŀ����
//		Day:		�趨��������Ϊ��֤����ͨ��ȡ1��
//		branch:		������0�򽵽���1
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
//		rvt[6]		tf����λ���ٶȣ�km��km/s��ͨ������Ŀ����л��㣬�۲�߶Ȳ���
void single_imp3(const double m0, const double t0, const double* rv0, int ID, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, double* rvt)
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
			tG = tG_estimate(ID, t0, coe0, Day, bra);//�ȹ���tG�ټ���Ŀ���ʵ��λ��
			flag = GetMov(t0 + tG, ID, phi, lambda);

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
			if (as < Re_km || as > Re_km + 1000.0)continue;

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
			if (am * (1.0 + em) > Re_km + 1000.0)continue;

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
	{
		////����ĩ��ʱ�̵�����״̬��ʵ�ʡ�����
		//for (int i = 0; i < 3; i++)
		//{
		//	rv0_[i] = rv0[i] * 1000.0;
		//	rv0_[i + 3] = (rv0[i + 3] + dv[i]) * 1000.0;
		//}
		//j2rv02rvf(rv0_, tf, rvf);
		//for (int i = 0; i < 6; i++)rvf[i] /= 1000, 0;

		//hf = sqrt(rvf[0] * rvf[0] + rvf[1] * rvf[1] + rvf[2] * rvf[2]) - Re_km;//ʵ�ʹ���߶ȣ���λ��km

		//if (hf > 500.0 || hf < 200.0)return;

		//rtar_F[0] = Re_km * cos(phi) * cos(lambda);//����Ŀ���λ������ת��
		//rtar_F[1] = Re_km * cos(phi) * sin(lambda);
		//rtar_F[2] = Re_km * sin(phi);
		//alphaGf = alphaG0_ + omegaE * tf;//tfʱ��GMST
		//rtar[0] = rtar_F[0] * cos(alphaGf) - rtar_F[1] * sin(alphaGf);
		//rtar[1] = rtar_F[0] * sin(alphaGf) + rtar_F[1] * cos(alphaGf);
		//rtar[2] = rtar_F[2];
		//for (int i = 0; i < 3; i++)
		//{
		//	rvt[i] = rtar[i] * (1.0 + hf / Re_km);
		//	rvt[i + 3] = rvf[i + 3];
		//}

		flag0 = 1;

	}


};


//Fix�ع�ϵkm��Geo����γ��rad
void Geo2Fix(double* Fix, const double* Geo)
{
	Fix[0] = Re_km * cos(Geo[1]) * cos(Geo[0]);
	Fix[1] = Re_km * cos(Geo[1]) * sin(Geo[0]);
	Fix[2] = Re_km * sin(Geo[1]);
}

//Eci����ϵkm��Fix�ع�ϵkm��tʱ��s
void Fix2Eci(double* Eci, const double* Fix, double t)
{
	double alpha_G = alpha_G0 + omegaE * t;
	Eci[0] = cos(alpha_G) * Fix[0] - sin(alpha_G) * Fix[1];
	Eci[1] = sin(alpha_G) * Fix[0] + cos(alpha_G) * Fix[1];
	Eci[2] = Fix[2];
}

//Eci����ϵkm��Geo����γ��rad��tʱ��s
void Geo2Eci(double* Eci, const double* Geo, double t)
{
	double Fix[3];
	Geo2Fix(Fix, Geo);
	Fix2Eci(Eci, Fix, t);
}

//�����ӳ���gamma����λrad
int cal_gamma(double& gamma, const double* robs, int tar_ID, double t)
{
	double rtar[3], Geo[2], dr[3];

	gamma = DPI / 2;
	if (V_Norm2(&robs[0], 3) < Re_km + 200.0)return 0;
	//if (V_Norm2(&robs[0], 3) > Re_km + 500.0)return 0;

	if (tar_ID < 201)
	{
		//Geo[0] = points[tar_ID - 1].Pos.longitude;
		//Geo[1] = points[tar_ID - 1].Pos.latitude;
	}
	else
	{
		int flag = GetMov(t, tar_ID, Geo[1], Geo[0]);
		if (flag == 0)return 0;
	}

	Geo2Eci(rtar, Geo, t);

	if (V_Dot(robs, rtar, 3) > 0.0)
	{
		for (int i = 0; i < 3; i++)dr[i] = robs[i] - rtar[i];
		double C = V_Dot(robs, dr, 3) / V_Norm2(robs, 3) / V_Norm2(dr, 3);
		if (C > 1)C = 1;
		gamma = acos(C);
		if (gamma <= 2.5 * D2R)return 1;
	}


	return 0;
}

//���ת��ʱ������������С��Ŀ��Ϊ�۲��
int shootfun1(int n, const double* x, double* fvec, int iflag, const double* para)
{
	double t0, tf, mdv, rv0[6], rvt[6], dv[3], mass;
	int flag, ID;
	tf = x[0];
	mdv = x[1];
	for (int i = 0; i < 6; i++)rv0[i] = para[i];
	mass = para[6] * exp(-fabs(mdv) * 1000.0 / g0 / Isp);
	t0 = para[7];
	ID = (int)para[8];

	for (int i = 0; i < 3; i++)dv[i] = mdv * rv0[i + 3] / V_Norm2(&rv0[3], 3);
	for (int i = 0; i < 3; i++)rv0[i + 3] += dv[i];

	PropagateODE45(mass, true, rv0, t0, t0 + tf, rvt);
	flag = cal_gamma(fvec[0], rvt, ID, t0 + tf);
	fvec[1] = fvec[0];

	return 1;
}


//����̫���߶Ƚǣ��Ƕȵ�λrad��ʱ���Ǵ�alpha_G0����ǰʱ�̵�����
int cal_sunangle(double& sunangle, double t, int tar_ID)
{
	const double JD0 = 2458849.5;//UTC 2020��1��1��0ʱ
	double alpha, alphaG, delta_sun, alpha_sun, JD, TJD, lam_M_sun, M_sun, lam_ecl, epsilon;
	double lambda, phi;

	if (tar_ID < 201)//�����Ǿ�ֹĿ��
	{
		//lambda = point_all[tar_ID - 1]->longitude;
		//phi = point_all[tar_ID - 1]->latitude;
	}
	else
	{
		int flag = GetMov(t, tar_ID, phi, lambda);
		if (flag == 0)return 0;
	}

	alphaG = alpha_G0 + omegaE * t;
	alpha = lambda + alphaG;

	JD = JD0 + t / 86400.0;
	TJD = (JD - 2451545.0) / 36525.0;

	lam_M_sun = 280.46 + 36000.771 * TJD;
	M_sun = DPI / 180.0 * (357.5277233 + 35999.05034 * TJD);
	lam_ecl = DPI / 180.0 * (lam_M_sun + 1.914666471 * sin(M_sun) + 0.019994643 * sin(2 * M_sun));
	epsilon = DPI / 180.0 * (23.439291 - 0.0130042 * TJD);

	delta_sun = asin(sin(epsilon) * sin(lam_ecl));
	alpha_sun = atan2(cos(epsilon) * sin(lam_ecl) / cos(delta_sun), cos(lam_ecl) / cos(delta_sun));

	sunangle = asin(sin(phi) * sin(delta_sun) + cos(phi) * cos(delta_sun) * cos(alpha_sun - alpha));
	if (sunangle < 20.0 * D2R)
		return 0;

	return 1;
}