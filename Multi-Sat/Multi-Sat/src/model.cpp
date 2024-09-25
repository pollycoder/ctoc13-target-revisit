#include "model.h"

//������Ƭ����R��ת��ʱ������ΪT�������ܵ��ٶ�����dv_sequence
//��Ƭ����R�а���n����Ƭ��˳��Ϊ��Ƭ�������˳��
//ת��ʱ������T�а���2(n-1)+1=2n-1��ʱ�̣���2i+1���͵�2i+2��Ԫ�طֱ�Ϊ����ƬR[i]ת����R[i+1]��ת����ʼʱ��ͽ���ʱ��(i=0,...,n-2)
//�����ʱ�䵥λΪday��������ٶȵ�λΪm/s

double estimate_T(const std::vector<int>& R, std::vector<double>& T, double initT)
{
	double dv = 0.0;
	T.clear();
	double tnow = initT;
	T.resize(2 * R.size() - 1);
	T[0]=(initT);
	for (int i= 1; i<R.size(); i++)
	{
		double Ts, Tf;
		dv += estimate_dv(R[i - 1], R[i], tnow, Ts, Tf, R.size()-i);
		T[2*i-1] = (Ts);
		T[2 * i]=(Tf);
		tnow = Tf;
	}
	return dv;
}

double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf)
{
	double ts = (tnow  + 5.0) * Day2Second;						//ת����ʼʱ��ts����Ϊ��ǰʱ��tnow��5���
	double tf_low = ts + 3.0 * 86400.0;							//tf���½磬ת����ʼʱ��ts���1��
	double tf_high = (tnow + 100.0)* Day2Second;					//tf���Ͻ磬��ǰʱ��tnow���30��
	
	//����ת����unvisited����Ƭ���ٶ���������ֵEstimatedv
	double a_s = debris_data[debris_now][2];
	double i_s = debris_data[debris_now][4];
	double Omega_s0 = debris_data[debris_now][5];
	double theta_s0 = debris_data[debris_now][6] + debris_data[debris_now][7];
	double a_f = debris_data[target][2];
	double i_f = debris_data[target][4];
	double Omega_f0 = debris_data[target][5];
	double theta_f0 = debris_data[target][6] + debris_data[target][7];
	
	double dOmega_s = -3.0 / 2.0 * J2 * (req / a_s) * (req / a_s) * sqrt(mu / a_s / a_s / a_s) * cos(i_s);
	double dOmega_f = -3.0 / 2.0 * J2 * (req / a_f) * (req / a_f) * sqrt(mu / a_f / a_f / a_f) * cos(i_f);

	//delta_Omega_low��delta_Omega_high��ʾ������ྭ��ı仯��Χ
	double delta_Omega_low  = Omega_f0 + dOmega_f * tf_low  - Omega_s0 - dOmega_s * tf_low;		//tf_low ʱ�̣���Ƭdebris_now����Ƭunvisited[j]������ྭ�Ĳ�ֵ
	double delta_Omega_high = Omega_f0 + dOmega_f * tf_high - Omega_s0 - dOmega_s * tf_high;	//tf_highʱ�̣���Ƭdebris_now����Ƭunvisited[j]������ྭ�Ĳ�ֵ

	int delta_Omega_low_mod2pi = floor(delta_Omega_low / 2.0 / pi);
	int delta_Omega_high_mod2pi = floor(delta_Omega_high / 2.0 / pi);
	double tf_temp;
	if (delta_Omega_low_mod2pi != delta_Omega_high_mod2pi)
		//��delta_Omega_low��delta_Omega_high֮�����2pi������������ѡ��tf_tempʹ��������ྭ��Ϊ2pi��������
		tf_temp = (std::max(delta_Omega_low_mod2pi, delta_Omega_high_mod2pi) * 2.0 * pi
			- Omega_f0 + Omega_s0) / (dOmega_f - dOmega_s);
	else if (std::min(delta_Omega_low - delta_Omega_low_mod2pi * 2.0 * pi, delta_Omega_low_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_low)
		< std::min(delta_Omega_high - delta_Omega_high_mod2pi * 2.0 * pi, delta_Omega_high_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_high))
		//��delta_Omega_high��ȣ�delta_Omega_low������2pi��������
		tf_temp = tf_low;
	else
		//��delta_Omega_low��ȣ�delta_Omega_high������2pi��������
		tf_temp = tf_high;

	Ts = ts / Day2Second;
	Tf = tf_temp / Day2Second;
	return Dv_ij(debris_now, target, Ts, Tf);
}

//NΪ�����û�вι�
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, int N)
{
	double ts = (tnow + 5.0) * Day2Second;						//ת����ʼʱ��ts����Ϊ��ǰʱ��tnow��5���
	double tf_low = ts + 3.0 * 86400.0;							//tf���½磬ת����ʼʱ��ts���1��

	double t_gap = (2952.0 - tnow) / N;
	if (t_gap > 115.0) t_gap = 115.0;
	
	double tf_high = (tnow + t_gap) * Day2Second;					//tf���Ͻ磬��ǰʱ��tnow���30��

	//����ת����unvisited����Ƭ���ٶ���������ֵEstimatedv
	double a_s = debris_data[debris_now][2];
	double i_s = debris_data[debris_now][4];
	double Omega_s0 = debris_data[debris_now][5];
	double theta_s0 = debris_data[debris_now][6] + debris_data[debris_now][7];
	double a_f = debris_data[target][2];
	double i_f = debris_data[target][4];
	double Omega_f0 = debris_data[target][5];
	double theta_f0 = debris_data[target][6] + debris_data[target][7];

	double dOmega_s = -3.0 / 2.0 * J2 * (req / a_s) * (req / a_s) * sqrt(mu / a_s / a_s / a_s) * cos(i_s);
	double dOmega_f = -3.0 / 2.0 * J2 * (req / a_f) * (req / a_f) * sqrt(mu / a_f / a_f / a_f) * cos(i_f);

	//delta_Omega_low��delta_Omega_high��ʾ������ྭ��ı仯��Χ
	double delta_Omega_low = Omega_f0 + dOmega_f * tf_low - Omega_s0 - dOmega_s * tf_low;		//tf_low ʱ�̣���Ƭdebris_now����Ƭunvisited[j]������ྭ�Ĳ�ֵ
	double delta_Omega_high = Omega_f0 + dOmega_f * tf_high - Omega_s0 - dOmega_s * tf_high;	//tf_highʱ�̣���Ƭdebris_now����Ƭunvisited[j]������ྭ�Ĳ�ֵ

	int delta_Omega_low_mod2pi = floor(delta_Omega_low / 2.0 / pi);
	int delta_Omega_high_mod2pi = floor(delta_Omega_high / 2.0 / pi);
	double tf_temp;
	if (delta_Omega_low_mod2pi != delta_Omega_high_mod2pi)
		//��delta_Omega_low��delta_Omega_high֮�����2pi������������ѡ��tf_tempʹ��������ྭ��Ϊ2pi��������
		tf_temp = (std::max(delta_Omega_low_mod2pi, delta_Omega_high_mod2pi) * 2.0 * pi
			- Omega_f0 + Omega_s0) / (dOmega_f - dOmega_s);
	//else if (std::min(delta_Omega_low - delta_Omega_low_mod2pi * 2.0 * pi, delta_Omega_low_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_low)
	//	< std::min(delta_Omega_high - delta_Omega_high_mod2pi * 2.0 * pi, delta_Omega_high_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_high))
	//	//��delta_Omega_high��ȣ�delta_Omega_low������2pi��������
	//	tf_temp = tf_low;
	else
		//��delta_Omega_low��ȣ�delta_Omega_high������2pi��������
		tf_temp = tf_high;

	Ts = ts / Day2Second;
	Tf = tf_temp / Day2Second;
	return Dv_ij(debris_now, target, Ts, Tf);
}

//�������Ƭi����Ƭj��ת����ʼʱ��Ϊts,ת�ƽ���ʱ��Ϊtf������ٶ�����
double Dv_ij(int i, int j, double ts, double tf)
{
	const double JD2S = 86400.0;													//��ת��
	const double muEarth = mu; // 3.98600441800e+14;								//m^3/s^2������������
	const double Re = req;// 6378137.0;												//m �������뾶

	double mjd0 = ts * Day2Second;
	double mjdf = tf * Day2Second;
	//���ݼ�debris_data
	//0��Ƭ��� 1* 2����볤��(m) 3ƫ���� 4���(rad) 
	//5������ྭ(rad) 6���ص����(rad) 7������(rad) 8�״ﷴ�����
	double a0 = debris_data[i][2];
	double af = debris_data[j][2];
	double e0 = debris_data[i][3];
	double ef = debris_data[j][3];
	double inc0 = debris_data[i][4];
	double incf = debris_data[j][4];
	double dOmega_s = -3.0 / 2.0 * J2 * (req / a0) * (req / a0) * sqrt(mu / a0 / a0 / a0) * cos(inc0);//��Ƭi��������ྭ�仯��
	double dOmega_f = -3.0 / 2.0 * J2 * (req / af) * (req / af) * sqrt(mu / af / af / af) * cos(incf);//��Ƭj��������ྭ�仯��
	double Omega0 = debris_data[i][5] + dOmega_s * mjdf;  //tfʱ�̵�ֵ��0�ǳ�����Ƭ��f�ǵ�����Ƭ
	double Omegaf = debris_data[j][5] + dOmega_f * mjdf;

	double omega0 = debris_data[i][6] + domega(i) * mjdf;
	double omegaf = debris_data[j][6] + domega(j) * mjdf;
	double dt = (mjdf - mjd0);

	double a = (a0 + af) / 2;
	double inc = (inc0 + incf) / 2;
	double V = sqrt(muEarth / a);

	double delta_Omega = fmod(Omegaf - Omega0, 2 * pi);
	if (delta_Omega > pi)
		delta_Omega -= 2 * pi;
	else if (delta_Omega <= -pi)
		delta_Omega += 2 * pi;

	double x = delta_Omega * sin(inc) * V;
	double y = (af - a0) / 2 / a * V;
	double z = (incf - inc0) * V;

	double Re_p0 = Re / a0 / (1 - e0 * e0);
	double dOmega0 = -1.5 * J2 * Re_p0 * Re_p0 * sqrt(muEarth / a0 / a0 / a0) * cos(inc0);
	double Re_pf = Re / af / (1 - ef * ef);
	double dOmegaf = -1.5 * J2 * Re_pf * Re_pf * sqrt(muEarth / af / af / af) * cos(incf);
	double dOmega = (dOmega0 + dOmegaf) / 2;

	double m = 7 * dOmega * sin(inc) * dt;
	double n = dOmega * tan(inc) * sin(inc) * dt;

	double sx = (2 * x + m * y + n * z) / (4 + m * m + n * n);
	double sy = (2 * m * x - (4 + n * n) * y + m * n * z) / (8 + 2 * m * m + 2 * n * n);
	double sz = (2 * n * x + m * n * y - (4 + m * m) * z) / (8 + 2 * m * m + 2 * n * n);

	double delta_x = m * sy + n * sz;
	double dV1, dV2;

	// ƫ��������
	double delta_ex = ef * cos(omegaf) - e0 * cos(omega0);
	double delta_ey = ef * sin(omegaf) - e0 * sin(omega0);
	double dVe = 0.5 * V * sqrt(delta_ex * delta_ex + delta_ey * delta_ey);
	dVe /= 2;
	dV1 = sqrt(sx * sx + sy * sy + sz * sz + dVe * dVe);
	dV2 = sqrt((x - sx - delta_x) * (x - sx - delta_x) + (y + sy) * (y + sy) + (z + sz) * (z + sz) + dVe * dVe);

	double dV;
	dV = dV1 + dV2;
	return dV;

}

double domega(int i)
{
	const double Re = req;// 6378137.0;													//m �������뾶
//���ݼ�debris_data
//0��Ƭ��� 1* 2����볤��(m) 3ƫ���� 4���(rad) 
//5������ྭ(rad) 6���ص����(rad) 7������(rad) 8�״ﷴ�����
	double a = debris_data[i][2];
	double e = debris_data[i][3];
	double inc = debris_data[i][4];
	double OMEGA = debris_data[i][5];
	double omega = debris_data[i][6];

	double p = a * (1 - e * e);
	double c2 = (Re / p) * (Re / p);
	double ci = cos(inc);
	double n = sqrt(mu / (a * a * a));
	double dOmega = -1.5 * J2 * c2 * n * ci;
	double dw = 0.75 * J2 * c2 * n * (5 * ci * ci - 1);
	return dw;
}

double Dv_All(double* T, int* R, int n)
{
	double dv = 0.0;
	for (int i = 0; i < n - 1; i++)
	{
		dv = dv + Dv_ij(R[i], R[i + 1], T[2 * i + 1], T[2 * i + 2]);
	}
	return dv;
}
double Dv_All(const std::vector<double>& T, int* R, int n)
{
	double dv = 0.0;
	for (int i = 0; i < n - 1; i++)
	{
		dv = dv + Dv_ij(R[i], R[i + 1], T[2 * i + 1], T[2 * i + 2]);
	}
	return dv;
}

double Dv_All(const std::vector<double>& T, const std::vector<int>& R, std::vector<double>& dv_sequence)
{

	int n = R.size();
	dv_sequence.resize(n - 1);
	for (int i = 0; i < n - 1; i++)
	{
		dv_sequence[i] = Dv_ij(R[i], R[i + 1], T[2 * i+1], T[2 * i + 2]);
	}
	return 0.0;
}