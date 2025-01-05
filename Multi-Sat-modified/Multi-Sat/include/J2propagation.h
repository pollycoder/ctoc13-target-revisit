// orbit_propagation.h

#ifndef ORBIT_PROPAGATION_H
#define ORBIT_PROPAGATION_H

#include "Constant.h"  // �����������壬�� D2R, Re_km ��
#include <cmath>
#include <iostream>
#include <cstring>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "OrbitFun.h"
#include "OrbitMath.h"


// ȫ�ֱ�������
extern double targets_latitude_longitute[20][2];

// ��������

// ���������꣨γ�ȣ����ȣ�ת��Ϊ�ع�ϵ����
void Geodetic2Fixed(double* Geodetic, double* Fixed);

// ���������꣨γ�ȣ����ȣ�ת��Ϊ���Ĺ���ϵ���꣬ʱ��Ϊ t
void Geodetic2J2000(double* Geodetic, double* R, double t);

// �������� 3x3 ����ĳ˻���C = A * B
void matrix_multi_only_three(double A[3][3], double B[3][3], double C[3][3]);

// ���� 3x3 ������ 3x1 �����ĳ˻���C = A * B
void martri_multiply_vector_only_3(double A[3][3], double* B, double* C);

// ����� J2000 ����ϵ���ع�ϵ����ת����ʱ��Ϊ t
void get_J20002FIX_Rotation(double t, double rotation_matrix[3][3]);

// �� 3x3 �������ת�ã��͵��޸ģ�
void get_trans(double rotation_matrix[3][3]);

// ����ӵع�ϵ�� J2000 ����ϵ����ת����ʱ��Ϊ t
void get_FIX2J2000_Rotation(double t, double rotation_matrix[3][3]);


// ���� J2 �Ŷ��Ĺ������ѧ���̣�λ���ٶ�����
int dyn_j2_rv(const double* x, double* dx);

// �߾��ȹ������ѧģ�ͣ�t Ϊ����ѧʱ�䣨�룩��x Ϊ״̬����
int Dynamic_model(double t, const double* x, double* dx, const double* hpara);

// ���� J2 �Ŷ��Ĺ�����֣��� t_start ���ֵ� t_end
int propagate_j2(double rv0[6], double rvf[6], double t_start, double t_end, double abstol_value = 1.0e-10, double restol_value = 1.0e-10);
//void propagate_j2linear(double* RV0, double* RVf, const double& t0, const double& tf);

void propagate_linearJ2(const  double* rv0, double* rvf, double t0, double tf);

// ��ȡ����Ŀ���ź�ʱ�� t �µĵ������꣨γ�ȣ����ȣ�
void get_target_geogetic(int id, double t, double* Geodetic);

// ��ȡĿ���� J2000 ����ϵ�µ�λ������
void get_target_R(int id, double t, double* RV);

// ���Թ�����ֺ���
int test_orbit_propgation();

// ���Ե�������ת������
int test_geodetic();

// ����ʱ���
inline std::time_t parseTimestamp(const std::string& timestamp) {
    std::tm tm = {};
    std::istringstream ss(timestamp);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
    return std::mktime(&tm);
}

// ��ȡ����ʱ�����2035��9��26��12��00��00������
double get_t_second(const std::string& timestamp);

// ��ȡ����ʱ����µ�Ŀ��λ��������γ�ȣ����ȣ�
inline void get_target_geogetic(const int& id, const std::string& timestamp, double* Geodetic) {
    double t = get_t_second(timestamp);
    get_target_geogetic(id, t, Geodetic);
}

// ��ȡ����ʱ����µ�Ŀ����Ĺ���ϵ���꣨�ٶ�δ��ֵ�������ʼ����
inline void get_target_R(const int& id, const std::string& timestamp, double* RV) {
    double t = get_t_second(timestamp);
    get_target_R(id, t, RV);
}

#endif  // ORBIT_PROPAGATION_H
