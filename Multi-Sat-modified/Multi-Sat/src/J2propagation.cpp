#include "J2propagation.h"
#include "Constant.h"
#include "ODE45.h"
#include "OrbitMath.h"
#include <chrono>


double targets_latitude_longitute[20][2] =
{
	{23.701 * D2R, 120.5 * D2R},
	{36.908 * D2R, 127.879 * D2R},
	{40.197 * D2R, 126.361 * D2R},
	{56.718 * D2R, 38.243 * D2R},
	{49.409 * D2R, 28.066 * D2R},
	{18.442 * D2R, 42.819 * D2R},
	{15.505 * D2R, 49.77 * D2R},
	{9.984 * D2R, 49.514 * D2R},
	{-24.539 * D2R, 32.108 * D2R},
	{43.923 * D2R, 23.521 * D2R},
	{37.951 * D2R, 33.445 * D2R},
	{35.402 * D2R, -116.512 * D2R},
	{36.107 * D2R, -77.997 * D2R},
	{31.315 * D2R, -83.652 * D2R},
	{4.773 * D2R, -72.428 * D2R},
	{-49.807 * D2R, -70.047 * D2R},
	{23.282 * D2R, 105.846 * D2R},
	{28.182 * D2R, 94.039 * D2R},
	{28.224 * D2R, 78.13 * D2R},
	{46.963 * D2R, -67.55 * D2R}
};



//���ϵת�ع�ϵ
//�������꣨γ�ȣ����ȣ�
void Geodetic2Fixed(double* Geodetic, double* Fixed)
{
	double ee = 6.69437999014E-3;
	double singeodetic = sin(Geodetic[0]);
	double xx = sqrt(1.0 - ee * singeodetic * singeodetic);
	
	double Amp = Re_km / xx ;
	double cos_Geodetic_0 = cos(Geodetic[0]);
	Fixed[0] = Amp * cos_Geodetic_0 * cos(Geodetic[1]);
	Fixed[1] = Amp * cos_Geodetic_0 * sin(Geodetic[1]);
	Fixed[2] = (Amp - Re_km * ee / xx) * singeodetic;
	return;
}



// AXB =C
void matrix_multi_only_three(double A[3][3],double B[3][3], double C[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			C[i][j] = 0.0;
			for (int k = 0; k < 3; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return;
}

//AXB = C
void martri_multiply_vector_only_3(double A[3][3], double* B, double *C)
{
	for (int i = 0; i < 3; i++)
	{
		C[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			C[i] += A[i][j] * B[j];
		}
	}
	return;
}

void get_J20002FIX_Rotation(double t,double rotation_matrix[3][3])
{
	//const double alpha_G0 = 3.2310939479887431; //rad

	double GST = alpha_G0 + omegaE * t;

	//��z��ת��GST Ϊ R2
	double cosGST = cos(GST);
	double sinGST = sin(GST);
	double R2[3][3] = {
		{cosGST,sinGST,0.},
		{-sinGST,cosGST,0.},
		{0,0,1.0}
	};

	//����¶�ϵ����2035-09-26 12:00:00.000
	double R1[3][3] = {
		{0.9999623332927, -0.00795648703,   -0.0034568483},
		{0.00795661805,  0.9999682804575,   2.40846051e-5},
		{0.0034565374,   -5.15882829e-5,   0.9999940533582468}
	};

	//����R2XR1�����ظ�rotation_matrix
	matrix_multi_only_three(R2, R1, rotation_matrix);

	return;
}

//��������ת�����ת�ã�Ҳ���������
void get_trans(double rotation_matrix[3][3])
{
	double new_mat[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			new_mat[i][j] = rotation_matrix[j][i];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			rotation_matrix[i][j] = new_mat[i][j];

	return;
}

void get_FIX2J2000_Rotation(double t, double rotation_matrix[3][3])
{
	get_J20002FIX_Rotation(t, rotation_matrix);
	get_trans(rotation_matrix);
}


//���ϵת���Ĺ���ϵ
//�������꣨γ�ȣ����ȣ�
void Geodetic2J2000(double* Geodetic, double* R, double t)
{
	double RF[3];
	Geodetic2Fixed(Geodetic, RF);

	double rotation_matrix[3][3];
	get_FIX2J2000_Rotation(t, rotation_matrix);

	martri_multiply_vector_only_3(rotation_matrix, RF, R);

	return;
}

//λ���ٶ�������j2����ѧ����
int dyn_j2_rv(const double* x, double* dx)
{
	double mu = mu_km_s;
	double j2 = J2;
	double radius = Re_km;

	double ri = 1.0 / sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	double ri2 = ri * ri;
	double muj2radius2r2 = 1.5 * mu * j2 * radius * radius * ri2;
	double frz = mu + muj2radius2r2 * (1.0 - 5.0 * x[2] * x[2] * ri2);
	double fr = 2.0 * muj2radius2r2;
	double ri3 = ri2 * ri;

	dx[0] = x[3];
	dx[1] = x[4];
	dx[2] = x[5];
	dx[3] = -x[0] * ri3 * frz;
	dx[4] = -x[1] * ri3 * frz;
	dx[5] = -x[2] * ri3 * (frz + fr);

	return 1;
}

//�߾��ȹ�������Һ�����tΪ����ѧʱMJD����λs��xΪλ���ٶȣ���λm��m/s
int Dynamic_model(double t, const double* x, double* dx, const double* hpara)
{
	/*double* para = (double*) hpara;*/

	dx[0] = x[3];
	dx[1] = x[4];
	dx[2] = x[5];

	double t_real = t;

	double rotation_matrix[3][3];

	//������ת����
	get_J20002FIX_Rotation(t, rotation_matrix);

	double CartesianR[3] = { x[0], x[1], x[2] };
	double CartesianV[3] = { x[3], x[4], x[5] };
	double FixedR[3];
	double FixedV[3];
	martri_multiply_vector_only_3(rotation_matrix, CartesianR, FixedR);
	martri_multiply_vector_only_3(rotation_matrix, CartesianV, FixedV);

	//����������������
	double rv[6] = { FixedR[0],FixedR[1],FixedR[2],FixedV[0],FixedV[1],FixedV[2] }, drv[6];
	dyn_j2_rv(rv, drv);
	double Accel[3] = { drv[3], drv[4], drv[5] };

	//�����
	get_trans(rotation_matrix);

	//����J2��+����
	double SumAccel[3] = { 0.0, 0.0, 0.0 };
	martri_multiply_vector_only_3(rotation_matrix, Accel, SumAccel);
	double dv_j2[3] = { SumAccel[0], SumAccel[1], SumAccel[2] };

	dx[3] = SumAccel[0];
	dx[4] = SumAccel[1];
	dx[5] = SumAccel[2];

	return 1;
}

int propagate_j2(double rv0[6], double rvf[6],  double t_start, double t_end, double abstol_value , double restol_value )
{
	int n = 6;
	memcpy(rvf, rv0, n * sizeof(double));
	double* abstol = new double[n];
	double* newwork = new double[n * 10];
	for (int i = 0; i < n; i++)  abstol[i] = abstol_value;

	//auto start_t = std::chrono::steady_clock::now();
	double para[1];
	int num;
	int flag = ode45(Dynamic_model, rvf, para, t_start, t_end, n, num, newwork, abstol, restol_value,0,-1,20.0);

	//auto end_t = std::chrono::steady_clock::now();
	//double duration_t = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 1000.0;
	//std::cout << "Time is " << duration_t << "s ; which is " << duration_t / 3600.0 << " hours" << std::endl;

	//std::cout<<rvf[0]<<" "<<rvf[1]<<" "<<rvf[2]<<" "<<rvf[3]<<" "<<rvf[4]<<" "<<rvf[5]<<std::endl;

	delete[] abstol;
	delete[] newwork;

	return flag;
}

void get_target_geogetic(int id, double t, double* Geodetic)
{
	if(0<=id && id<20)
	{
		Geodetic[0] = targets_latitude_longitute[id][0];
		Geodetic[1] = targets_latitude_longitute[id][1];
	}
	else if (id == 20) // ship
	{
		// ��ʼ�ͽ�����γ�Ⱥ;��ȣ���λ���ȣ�
		double lat_start_deg = 10.642;
		double lon_start_deg = 53.558;
		double lat_end_deg = 6.582;
		double lon_end_deg = 77.963;

		// ����ʼ�ͽ�����γ�Ⱥ;���ת��Ϊ����
		double lat1 = lat_start_deg * D2R;
		double lon1 = lon_start_deg * D2R;
		double lat2 = lat_end_deg * D2R;
		double lon2 = lon_end_deg * D2R;

		// ��ʼ�ͽ�����ʱ�䣨��λ���룩
		double t_start = 0.0;
		double t_end = 86400.0 * 2.0; // 86400�� * 2��

		// ȷ�� t ����Ч��Χ��
		if (t < t_start) t = t_start;
		if (t > t_end) t = t_end;

		// ����ʱ�����
		double ratio = (t - t_start) / (t_end - t_start);

		// ���������յ�֮��ĽǾ��� delta_sigma
		double delta_lon = lon2 - lon1;
		double sin_lat1 = sin(lat1);
		double cos_lat1 = cos(lat1);
		double sin_lat2 = sin(lat2);
		double cos_lat2 = cos(lat2);
		double sin_delta_lon = sin(delta_lon);
		double cos_delta_lon = cos(delta_lon);

		double delta_sigma = acos(sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lon);

		// �����ʼ��λ�� alpha1
		double y = cos_lat2 * sin_delta_lon;
		double x = cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_delta_lon;
		double alpha1 = atan2(y, x);

		// ���㵱ǰ�Ǿ��� delta_sigma_f
		double delta_sigma_f = ratio * delta_sigma;

		// ���㵱ǰλ�õ�γ�� phi_f
		double sin_delta_sigma_f = sin(delta_sigma_f);
		double cos_delta_sigma_f = cos(delta_sigma_f);

		double sin_phi_f = sin_lat1 * cos_delta_sigma_f + cos_lat1 * sin_delta_sigma_f * cos(alpha1);
		double phi_f = asin(sin_phi_f);

		// ���㵱ǰλ�õľ��� lambda_f
		y = sin(alpha1) * sin_delta_sigma_f * cos_lat1;
		x = cos_delta_sigma_f - sin_lat1 * sin_phi_f;
		double lambda_f = lon1 + atan2(y, x);

		// �����ȹ淶���� [-��, ��]
		lambda_f = fmod(lambda_f + 3 * DPI, 2 * DPI) - DPI;

		// ���ؼ�����
		Geodetic[0] = phi_f;    // γ�ȣ����ȣ�
		Geodetic[1] = lambda_f; // ���ȣ����ȣ�
	}
	else
	{
		std::cout << "wrong number!" << std::endl;
	}

	return;
}


void get_target_R(int id, double t, double* RV)
{

	double Geodetic[2];
	get_target_geogetic(id, t, Geodetic);
	Geodetic2J2000(Geodetic, RV, t);
	//std::cout << "t" << t << "  ;" << RV[0] << " " << RV[1] << " " << RV[2] << std::endl;
}

//���Թ������
int test_orbit_propgation()
{
	
	double rv0[6] = { -1197.4673034427, 6791.17454760102 ,0.0 ,3.85623032827321,0.679957450700828,6.51685061874543 };
	double rvf[6];
	double t_start = 0.0;
	double t_end = 86400.0*2.0;

	double abstol_value = 1.0e-10;
	double restol_value = abstol_value;

	int flag = propagate_j2(rv0, rvf, t_start, t_end, abstol_value, restol_value);

	std::cout<<rvf[0]<<" "<<rvf[1]<<" "<<rvf[2]<<" "<<rvf[3]<<" "<<rvf[4]<<" "<<rvf[5]<<std::endl;

	return flag;
}

int test_geodetic()
{
	double lambda0 = 23.701 * D2R; //Ŀ���ĵ���γ�ȣ�degree
	double phi0 = 120.5 * D2R; //Ŀ���ĵ��澭�ȣ�degree

	double Geodetic[2] = { lambda0 , phi0 };

	double RV[3];
	double t = 0.0;
	Geodetic2J2000(Geodetic, RV, t);

	std::cout <<"t" << t << "  ;" << RV[0] << " " << RV[1] << " " << RV[2] << std::endl;


	t = 86400.0*2.0;

	Geodetic2J2000(Geodetic, RV, t);
	std::cout << "t" << t << "  ;" << RV[0] << " " << RV[1] << " " << RV[2] << std::endl;

	t = 86400.0 * 1.0;
	get_target_R(20, t, RV);

	std::cout << "t: " << t << "  ;" << RV[0] << " " << RV[1] << " " << RV[2] << std::endl;


	return 1;
}


double get_t_second(const std::string& timestamp) {
	// ����Ŀ��ʱ��
	std::tm targetTime = {};
	targetTime.tm_year = 2035 - 1900; // ��ݴ�1900��ʼ����
	targetTime.tm_mon = 8; // �·ݴ�0��ʼ������9����8
	targetTime.tm_mday = 26;
	targetTime.tm_hour = 12;
	targetTime.tm_min = 0;
	targetTime.tm_sec = 0;

	// ��Ŀ��ʱ��ת��Ϊʱ���
	std::time_t targetTimestamp = std::mktime(&targetTime);

	// ��������ʱ���
	std::time_t inputTimestamp = parseTimestamp(timestamp);

	// ʹ��chrono������ֵ
	auto inputTimePoint = std::chrono::system_clock::from_time_t(inputTimestamp);
	auto targetTimePoint = std::chrono::system_clock::from_time_t(targetTimestamp);

	auto duration = std::chrono::duration_cast<std::chrono::seconds>(inputTimePoint - targetTimePoint);

	double result = static_cast<double>(duration.count());

	return result;
}