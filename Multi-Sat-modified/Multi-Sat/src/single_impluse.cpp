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

	Nmax = floor(10800.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));//���ת��Ȧ��(3h)
	Nmin = 0;//floor((double(Day) - 1.0) * 86400.0 / (2.0 * DPI * sqrt(a0 * a0 * a0 / mu_km_s)));// +1;//����ת��Ȧ��

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
//TODO��ȥ������branch��NR�ı�����Ҳ��Ҫ�������������ɹ�������չ
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
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch) 
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

	rv2coe(flag0, coe0, rv0, mu_km_s);//�������ת���ɹ������

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

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	//for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//�����ͽ�����֧
		int bra = branch;
		NRi = NR;
		{
			// ZZ �������
			if (phi > inc0 + 0.1) {
				flag0 = 0;
				return;
			}
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
			
			//��ʱam���������i,Omega���䣬���������ƽ������������ƽ⣩
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
			double temp1 = mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0));
			double temp2 = mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00));
			if (temp1 < 0 || temp2 < 0) {
				flag0 = 0;
				return;
			}
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));
			dv_min = dv_mod;

			// �ж��Ƿ�Ϊnan
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

	omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

	c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
	c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
	alphaGt = Omega0 + atan2(c1, c2);
	while (alphaGt < alphaG0_)alphaGt += D2PI;
	while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
	tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / omegaE;
}


//Tg �Ǿ����೤ʱ�䣬����ʱ��
void modify_tg(const double t0, const double* rv0, double phi, double& tG)
{
	//�ٴ�����ʱ��
	double rv_temp[6];
	memcpy(rv_temp, rv0, 6 * sizeof(double));
	propagate_linearJ2(rv0, rv_temp, t0, tG + t0);

	int flag_coe;
	double coe[6];
	rv2coe(flag_coe, coe, rv_temp, mu * 1e-9);  //���������

	double u0 = coe[4] + coe[5]; u0 = NiceAngle(u0);//���//��Χ��0 ��2pi
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
	//��Χ�������pi 
	double du1 = NiceAngle(u1 - u0);
	double du2 = NiceAngle(u2 - u0);
	if (du1 > DPI) du1 = du1 - 2 * DPI;
	if (du2 > DPI) du2 = du2 - 2 * DPI;

	//startʱ�̵�λ���ٶ�
	double r_start_staellite[3] = { rv_temp[0],rv_temp[1],rv_temp[2] };
	double v_start_staellite[3] = { rv_temp[3],rv_temp[4],rv_temp[5] };

	//Ѱ�����µ�γ�ȵ�һ����Ŀ���غϵ�����ʱ��
	double _omega_ = V_Norm2(v_start_staellite, 3) / V_Norm2(r_start_staellite, 3);
	double t1 = du1 / _omega_;
	double t2 = du2 / _omega_;
	double t = fabs(t1) < fabs(t2) ? t1 : t2;
	double  t_access = tG + t;

	if (t_access > 0)
    	tG = t_access;
}

//���������(ָ����������С�����:������Lambert���)
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
//		tG:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp_zzmodified(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tG, double* dv, const int& branch)
{
	double lambda, phi;//Ŀ���ľ�γ�ȣ�ת���ɻ���
	double a0, e0, inc0, Omega0, omega0, f00, E00, M00, gamma00, coe0[6], ocoe0[6];//��ʼ�������
	double c1, c2;//�м����
	double alphaG0_, alphaGt;//GMST
	double omega_ft;
	int flag;

	flag0 = 0;

	rv2coe(flag0, coe0, rv0, mu_km_s);//�������ת���ɹ������

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

	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]


	int bra = branch;

	double inc_lower90 = inc0;
	if (inc_lower90 > DPI / 2.0) inc_lower90 = DPI - inc_lower90;

	if(phi > inc0 + 0.06)
	{
		flag0 = 0;
		return;
	}

	//һ���Ž�������
	if (phi > inc0 -0.06) {
		//ZZ TODO: ��Ҫ����86400 ��������

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

	 //���㾫ȷ����ʱ�̵�R����
	double Geodetic[2] = { phi0,lambda0 };
	double RV_Target[6];
	Geodetic2J2000(Geodetic, RV_Target, tG + t0);
	V_Multi(RV_Target, RV_Target, a0 / V_Norm2(RV_Target, 3), 3);

	//////////////////////////////////////����Lambert����////////////////////////
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
	am = (V_Norm2(rv0, 3) + V_Norm2(RV_Target, 3) + V_Norm2(vt0, 3)) / 4.0; //��С��Բת�ƹ���볤��
	Nmax = (int)floor(tG / (6.283185307179586476925286766559 * sqrt(am * am * am / mu_km_s))) + 1; //�����ܵ�ת��Ȧ��
	int flag_total = 0;

	double dt_min = 0.0;

	for (int n = 0; n <= Nmax; n++) //Ȧ����0������������٣����������ٶ���С�Ľ��
	{
		for (int j = 0; j < 2; j++) //���ڶ�Ȧ���⣬�����ҷ�֦������
		{
			double dt = 0.0;
			//for(double dt = -2000.0; dt < 2010.0 ; dt+=100.0)
			{
				
				//���㾫ȷ����ʱ�̵�R����
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

							if (flag0 != 1 || e >= 1.0) continue;//�ų�˫����������
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

	tG = tG + dt_min; //�˴��ܹؼ������ܱ��tG + dt_min������ͨ���Ż�����
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

//���������(ָ����������С�����:������Lambert���)
void lambert_dv(int& flag, const double t0, const double* rv0,  double* dv,  double tf, int id)
{
	flag = 0;
	if (tf - t0 < 10.0) return;

	double coe0[6];
	int flag_orbit;
	rv2coe(flag_orbit, coe0, rv0, mu_km_s);//�������ת���ɹ������
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
	am = a0; //��С��Բת�ƹ���볤��
	Nmax = (int)floor((tf-t0) / (6.283185307179586476925286766559 * sqrt(am * am * am / mu_km_s))) + 1; //�����ܵ�ת��Ȧ��


	for (int n = 0; n <= Nmax; n++) //Ȧ����0������������٣����������ٶ���С�Ľ��
	{
		for (int j = 0; j < 2; j++) //���ڶ�Ȧ���⣬�����ҷ�֦������
		{
			
			for(int way = 0; way < 2; way++)
			{
				Geodetic2J2000(Geodetic, RV_Target, tf );
				V_Multi(RV_Target, RV_Target, a0 / V_Norm2(RV_Target, 3), 3);
				int flag0;
				lambert(v0, v1, a, e, rv0, RV_Target, tf  - t0, unith, flag0, mu_km_s, way, n, j, 200);

				if (flag0 != 1 || e >= 1.0) continue;//�ų�˫����������
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

	double lon1 = 53.558 * D2R;
	double lon2 = 77.963 * D2R;
	double omega_ship = (lon2 - lon1) / 2.0 / 86400.0;


	flag0 = 0;

	rv2coe(flag0, coe0, rv0, mu_km_s);//�������ת���ɹ������

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

	k1 = (1.0 - e0 * e0) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) * (sin(f00 - gamma00) + sin(f00) * cos(gamma00) / (1.0 + e0 * cos(f00)));
	k2 = (1.0 - e0 * e0) / 2.0 * (cos(f00 - gamma00) + cos(E00) * cos(gamma00)) / (e0 * cos(f00 - gamma00) + cos(gamma00));
	k3 = sqrt((1.0 - e0 * e0) * (1.0 - e0 * e0) * (1.0 - e0 * e0)) / (2.0 * e0) / (e0 * cos(f00 - gamma00) + cos(gamma00)) *
		(sin(f00 - gamma00) + (2.0 * e0 * sin(gamma00) + sin(f00) * cos(gamma00)) / (1.0 + e0 * cos(f00)));

	alphaG0_ = alpha_G0 + omegaE * t0;//��ʼʱ��GMST
	while (alphaG0_ > D2PI)alphaG0_ -= D2PI;//ע�ⷶΧ[0, 2PI]

	dv_min = 1.0e10;
	tmin = 0.0;
	//for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//�����ͽ�����֧
		int bra = branch;
		NRi = NR;
		{
			if (phi > inc0) {
				flag0 = 0;
				return;
			}
			omega_ft = (bra == 0) ? asin(sin(phi) / sin(inc0)) : -asin(sin(phi) / sin(inc0)) + DPI;//����omega+ft�������չ���Ĳ��������ǹ����û��

			c1 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (cos(omega_ft) * sin(lambda) - sin(omega_ft) * cos(inc0) * cos(lambda));
			c2 = -1 / sqrt(1 - sin(phi) * sin(phi)) * (-sin(omega_ft) * cos(inc0) * sin(lambda) - cos(omega_ft) * cos(lambda));
			alphaGt = Omega0 + atan2(c1, c2);
			while (alphaGt < alphaG0_)alphaGt += D2PI;
			while (alphaGt > alphaG0_ + D2PI)alphaGt -= D2PI;
			tG = ((alphaGt - alphaG0_) + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE + omega_ship);

			fs = (bra == 0) ? asin(sin(phi) / sin(inc0)) - omega0 : -asin(sin(phi) / sin(inc0)) - omega0 + DPI;

			//������(����as)�����η��̣�����fabs(as-a0)��С������
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


			//����J2ģ�����
			C_J2 = 1.5 * J2 * Re_km * Re_km * sqrt(mu_km_s) * pow(as, -3.5);
			omega_J2_ = C_J2 * (2.0 - 2.5 * sin(inc0) * sin(inc0)) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			Omega_J2_ = -C_J2 * cos(inc0) / (1.0 - e0 * e0) / (1.0 - e0 * e0);
			M_J2_ = C_J2 * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / pow((1.0 - e0 * e0), 1.5);
			tG_J2 = (alphaGt - alphaG0_ + 2.0 * DPI * (double(Day) - 1.0)) / (omegaE + omega_ship - Omega_J2_);

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

			//��ʱam���������i,Omega���䣬���������ƽ������������ƽ⣩
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

			//ƽ��ת��Ϊ˲��������ƽ˲��ת��������ʹ�õ�ReΪm��������Ҳ��ʱ�ĳ�m������
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

			//�����ʼʱ�̵���������
			p0 = a0 * (1 - e0 * e0);
			dv_mod = sqrt(mu_km_s * (2.0 * (1.0 + e0 * cos(f00)) / p0 - 1.0 / (ocoe[0] / 1000.0))) - sqrt(mu_km_s / p0 * (1.0 + e0 * e0 + 2.0 * e0 * cos(f00)));
			dv_min = dv_mod;

			// �ж��Ƿ�Ϊnan
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
	//for (NRi = Nmin; NRi <= Nmax; NRi++)
	{
		//for (int bra = 0; bra < 2; bra++)//�����ͽ�����֧
		int bra = branch;
		NRi = NR;
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
			//if (as < Re_km + 200.0 || as > Re_km + 1000.0)continue;

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
			if (am < Re_km + 200.0) {
				flag0 = 0;
				return;
			}

			//��ʱam���������i,Omega���䣬���������ƽ������������ƽ⣩
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

			if ((double)sign * dv_mod < 0.0) {
				flag0 = 0;
				return;//�������巽���ѡ��
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
//TODO: �޸��Ż�ʱ��propagate��ʹ��J2����ģ��
//�Ż����������㶯���
//		dv[3]�����壬ÿ�������Ŷ���30m/s
//		tf��ʵ�ʵ���ĩʱ�̣���6h��ʼ�Ŷ�������Ϊ[0,12h]
//�Ż�ָ�꣺
//		dv�������С
//������
//		t0����ʼʱ��
//		RV0[6]����ʼλ���ٶ�
//		tf���Ż�ǰ����ĩʱ�̣��̶�3h
//		dv[3]���Ż�ǰ������
//		target_id��Ҫ�۲��Ŀ����ţ�0-20��
//Լ����
//		��֤Ŀ��ɼ�
//		����߶�200-1000km
void perturbation(double* dv, double& tf, const std::vector<double>& X) {
	tf += (X[3] - 0.5) * 1.0e4;								
	dv[0] += (X[0] - 0.5) * 1.0e2;
	dv[1] += (X[1] - 0.5) * 1.0e2;
	dv[2] += (X[2] - 0.5) * 1.0e2;
}

//----------------------------------------------
// �������������� tf�������Ӧ�� angle
//----------------------------------------------
double calcAngle(int id,
	double tf,
	const double RV1[6],
	double t0)
{
	// 1. ����Ŀ�������
	double target_R[3];
	get_target_R(id, tf, target_R);

	// 2. ���� RV1 ������ tf �õ� RVf
	double RVf[6];
	propagate_linearJ2(RV1, RVf, t0, tf);

	// 3. �����ʱ�����߽Ƕ�
	//double angle = angle_target_visible(RVf, target_R, 0.0);
	double angle = angle_earth_view_angle(RVf, target_R);
	return angle;
}

//----------------------------------------------
// �ƽ�ָ��������ĺ������� [ta, tb] ������ angle ����Сֵ
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
	//�ȴ�ɸ
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

	// �ƽ�ָ��
	const double phi = (1.0 + std::sqrt(5.0)) / 2.0;

	// ��ʼ�����ڵ��������Ե�
	double c = tb - (tb - ta) / phi;
	double d = ta + (tb - ta) / phi;

	// �������ȴ����趨�ľ���ʱ��������С����
	while ((tb - ta) > tolerance)
	{
		double angleC = calcAngle(id, c, RV1, t0);
		double angleD = calcAngle(id, d, RV1, t0);

		// ����� c ���� d ���� angle ��˵�����Ž��� [c, tb] ��
		// ������ [ta, d] ��
		if (angleC > angleD)
		{
			// [c, tb] ����
			ta = c;
			c = d;
			d = ta + (tb - ta) / phi;
		}
		else
		{
			// [ta, d] ����
			tb = d;
			d = c;
			c = tb - (tb - ta) / phi;
		}
	}

	// ���������е���Ϊ�������Ž�
	return 0.5 * (ta + tb);
}


//----------------------------------------------
// ���������������� [ta, tb] ����ȫ������
// ���躯������������ֲ���Сֵ�������ȣ���
// ���ȴֲ������ٶԿ��������ûƽ�ָ��Ҿֲ���ֵ�����ѡȫ����С
//----------------------------------------------
double global_min_search_two_valleys(
	int id,
	const double RV1[6],
	double t0,
	double ta,
	double tb,
	double coarseStep = 10.0,   // �ֲ����������ɸ����������
	double localTol = 1.0e-1  // �ƽ�ָ��������ݲ�
)
{
	// 1) �ֲ���
	std::vector<double> xs;
	for (double x = ta; x <= tb; x += coarseStep)
	{
		xs.push_back(x);
	}
	if (xs.back() < tb)
	{
		xs.push_back(tb); // ��֤�� tb ������ȥ
	}

	// 2) ����ÿ��������ĺ���ֵ
	std::vector<double> fs(xs.size());
	for (size_t i = 0; i < xs.size(); ++i)
	{
		fs[i] = calcAngle(id, xs[i], RV1, t0);
	}

	// 3) ���ݲ���������ҿ��ܴ��ڡ��ֲ��ȡ�������
	//    ��� f[i-1] >= f[i] <= f[i+1]������ [xs[i-1], xs[i+1]] ������м�Сֵ
	//    ע��߽�������Ե�������
	std::vector<std::pair<double, double>> candidateIntervals;
	for (size_t i = 1; i + 1 < xs.size(); ++i)
	{
		double fPrev = fs[i - 1];
		double fCurr = fs[i];
		double fNext = fs[i + 1];

		// �ж� "�Ƚ�����" ���� ">= <= >="
		if (( fCurr <= fPrev && fCurr <= fNext) || i ==1 || i+1 == xs.size()-1)
		{
			double left = xs[i - 1];
			double right = xs[i + 1];
			candidateIntervals.push_back({ left, right });
		}
	}

	// ������������Ƚϴ󣬻��ߺ����Ĺ������ڱ�Ե��Ҳ����û��ʶ��
	// ����������� [ta, tb] Ҳ�ӽ�������
	if (candidateIntervals.empty())
	{
		candidateIntervals.push_back({ ta, tb });
	}

	// 4) �ûƽ�ָ���������Щ��ѡ�������Ҿֲ���С
	double bestX = ta;
	double bestF = std::numeric_limits<double>::infinity();

	for (auto& seg : candidateIntervals)
	{
		double segA = seg.first;
		double segB = seg.second;

		// ע�⣺�� segA==segB �� segA>segB(�������)Ҫ������
		if (segA > segB) std::swap(segA, segB);
		if (std::fabs(segB - segA) < 1e-12)
		{
			// ����̫С���������ֶ�����
			continue;
		}

		// �����������������������
		double xLocMin = golden_section_search(id, RV1, t0, segA, segB, localTol);
		double fLocMin = calcAngle(id, xLocMin, RV1, t0);

		// ����ȵ�ǰ���Ż�С�������
		if (fLocMin < bestF)
		{
			bestF = fLocMin;
			bestX = xLocMin;
		}
	}

	// ���� (x, f)
	return bestX;
}


double obj_func_shooting(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	double* para = static_cast<double*>(f_data);

	// ��ֵ
	const double t0 = para[0];
	double tf = para[1];
	double dv[3] = { para[2], para[3], para[4] };
	const double RV0[6] = { para[5], para[6], para[7], para[8], para[9], para[10] };
	const int id = static_cast<int>(para[11]);

	perturbation(dv, tf, X);

	double penalty_time = 0.0;
	if (tf < t0) {
		//std::cout << "tf < t0��������Ҫ��" << std::endl;
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


	// �ö�����ƣ�Ϊ�˸�J2��������������20km������
	double peri = a * (1 - e) - Re_km;
	double apo = a * (1 + e) - Re_km;
	double penalty_orbit = 0.0;
	penalty_orbit += fabs(std::min(peri, 200.0) - 200.0);
	penalty_orbit += fabs(std::max(apo, 1000.0) - 1000.0);


	//propagate_j2(RV1, RVf, t0, tf, 1e-10, 1e-10);
	//propagate_linearJ2(RV1, RVf, t0, tf);
	//rv02rvf(flag, RVf, RV1, tf - t0, mu_km_s);


	//ZZ Test
	// ���ûƽ�ָ��������õ�ʹ angle ��С��ʱ�� tf_opt
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

//// TODO���Ż�֮���ø߾���propagate�ٻ�һ�η�
//void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const int& target_id, const int& NR, const int& branch) {
//	flag = 0;
//	double target_geo[2];
//	double lambda, phi;  //lambda ���ȣ�phi γ��
//	if (target_id != 20) {
//		get_target_geogetic(target_id, tf, target_geo);
//		lambda = target_geo[1];
//		phi = target_geo[0];
//	}
//	else { // ��
//		get_target_geogetic(target_id, 2.0 * 86400.0, target_geo);
//		double phi_end = target_geo[0];
//		get_target_geogetic(target_id, 0.0, target_geo); //TODO�� ZZ ��ʱ�Ƿ���� Ȧ��Ӧ���Ǵӿ�ͷ�����ڵ�Ȧ����������ʱ�����Ƹ�����ɣ�
//		double phi_start = target_geo[0];
//
//		// γ�ȸ���ȷһЩ����һȦ105min�ƣ���γ�ȵȷ�
//		double dphi = (phi_start - phi_end) * double(NR) * 6300.0 / 172800.0;
//		phi = phi_start + dphi;
//		lambda = target_geo[1];
//	}
//	
//	double m0 = 1000.0, mf; 
//
//	//ע���Ÿ����ĵķ����tf�Ƿ���ʱ��������ʱ��
//	//���ڵ���Ŀ�꣬�ſ��о�γ�ȣ��ȼ�����ļнǷ�Χ���ٸ�����������Ⱥ�γ�ȶ������񣬴���ѡ������͵�
//	//���ڴ���γ���Ѿ��ȽϽ���
//	if (target_id != 20) {
//		// ������ĽǷ�Χbeta
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
//	else { // ��
//		// ������ĽǷ�Χbeta
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
//	nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 10000);		//�����
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


// TODO���Ż�֮���ø߾���propagate�ٻ�һ�η�
void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const int& target_id, const int& NR, const int& branch) {
	flag = 0;
	double target_geo[2];
	double lambda, phi;  //lambda ���ȣ�phi γ��
	if (target_id != 20) {
		get_target_geogetic(target_id, tf, target_geo);
		lambda = target_geo[1];
		phi = target_geo[0];
	}
	else { // ��
		get_target_geogetic(target_id, 2.0 * 86400.0, target_geo);
		double phi_end = target_geo[0];
		get_target_geogetic(target_id, 0.0, target_geo); //TODO�� ZZ ��ʱ�Ƿ���� Ȧ��Ӧ���Ǵӿ�ͷ�����ڵ�Ȧ����������ʱ�����Ƹ�����ɣ�
		double phi_start = target_geo[0];

		// γ�ȸ���ȷһЩ����һȦ105min�ƣ���γ�ȵȷ�
		double dphi = (phi_start - phi_end) * double(NR) * 6300.0 / 172800.0;
		phi = phi_start + dphi;
		lambda = target_geo[1];
	}

	double m0 = 1000.0, mf;

	//ע���Ÿ����ĵķ����tf�Ƿ���ʱ��������ʱ��
	//���ڵ���Ŀ�꣬�ſ��о�γ�ȣ��ȼ�����ļнǷ�Χ���ٸ�����������Ⱥ�γ�ȶ������񣬴���ѡ������͵�
	//���ڴ���γ���Ѿ��ȽϽ���
	if (target_id != 20) {
		// ������ĽǷ�Χbeta
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
	else { // ��
		// ������ĽǷ�Χbeta
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
	nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 100000);		//�����

	//PSO(obj_func_shooting, f_data, X, impulse, X.size(), 0, 100000, 2000000);		//�����



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
	double lambda, phi;  //lambda ���ȣ�phi γ��
	if (target_id != 20) {
		get_target_geogetic(target_id, t0, target_geo);
		lambda = target_geo[1];
		phi = target_geo[0];
	}
	else { // ��
		get_target_geogetic(target_id, t0, target_geo);
		lambda = target_geo[1];
		phi = target_geo[0];
	}

	//double m0 = 1000.0, mf;
	//double tg = 0.0; // ע���Ǿ�����ʱ�䣬����ʱ��
	//double dv[3];
	//ע��tf�Ƿ���ʱ��������ʱ��
	// {
	//	// ������ĽǷ�Χbeta
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
		nlopt_main(obj_func_shooting, f_data, X, impulse, X.size(), 0, 1000);		//�����
		perturbation(dv_info_temp.dv, dv_info_temp.tf, X);

		//��֤�Ƿ��ܿɼ�����ע��Ҫ��������
		double RV0_temp[6];
		memcpy(RV0_temp, RV0, 6 * sizeof(double));
		for (int i = 0; i < 3; i++) RV0_temp[i + 3] += dv_info_temp.dv[i];
		double tf_60 = round(dv_info_temp.tf / 60) * 60.0;  //ȡ�������ӣ��������趨һ��
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