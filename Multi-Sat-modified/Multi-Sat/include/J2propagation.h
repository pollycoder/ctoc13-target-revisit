// orbit_propagation.h

#ifndef ORBIT_PROPAGATION_H
#define ORBIT_PROPAGATION_H

#include "Constant.h"  // 包含常量定义，如 D2R, Re_km 等
#include <cmath>
#include <iostream>
#include <cstring>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "OrbitFun.h"
#include "OrbitMath.h"


// 全局变量声明
extern double targets_latitude_longitute[20][2];

// 函数声明

// 将地理坐标（纬度，经度）转换为地固系坐标
void Geodetic2Fixed(double* Geodetic, double* Fixed);

// 将地理坐标（纬度，经度）转换为地心惯性系坐标，时间为 t
void Geodetic2J2000(double* Geodetic, double* R, double t);

// 计算两个 3x3 矩阵的乘积：C = A * B
void matrix_multi_only_three(double A[3][3], double B[3][3], double C[3][3]);

// 计算 3x3 矩阵与 3x1 向量的乘积：C = A * B
void martri_multiply_vector_only_3(double A[3][3], double* B, double* C);

// 计算从 J2000 坐标系到地固系的旋转矩阵，时间为 t
void get_J20002FIX_Rotation(double t, double rotation_matrix[3][3]);

// 对 3x3 矩阵进行转置（就地修改）
void get_trans(double rotation_matrix[3][3]);

// 计算从地固系到 J2000 坐标系的旋转矩阵，时间为 t
void get_FIX2J2000_Rotation(double t, double rotation_matrix[3][3]);


// 考虑 J2 扰动的轨道动力学方程，位置速度描述
int dyn_j2_rv(const double* x, double* dx);

// 高精度轨道动力学模型，t 为动力学时间（秒），x 为状态向量
int Dynamic_model(double t, const double* x, double* dx, const double* hpara);

// 考虑 J2 扰动的轨道积分，从 t_start 积分到 t_end
int propagate_j2(double rv0[6], double rvf[6], double t_start, double t_end, double abstol_value = 1.0e-10, double restol_value = 1.0e-10);
//void propagate_j2linear(double* RV0, double* RVf, const double& t0, const double& tf);

void propagate_linearJ2(const  double* rv0, double* rvf, double t0, double tf);

// 获取给定目标编号和时间 t 下的地理坐标（纬度，经度）
void get_target_geogetic(int id, double t, double* Geodetic);

// 获取目标在 J2000 坐标系下的位置向量
void get_target_R(int id, double t, double* RV);

// 测试轨道积分函数
int test_orbit_propgation();

// 测试地理坐标转换函数
int test_geodetic();

// 解析时间戳
inline std::time_t parseTimestamp(const std::string& timestamp) {
    std::tm tm = {};
    std::istringstream ss(timestamp);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
    return std::mktime(&tm);
}

// 获取任意时间距离2035年9月26日12：00：00的秒数
double get_t_second(const std::string& timestamp);

// 获取任意时间戳下的目标位置向量（纬度，经度）
inline void get_target_geogetic(const int& id, const std::string& timestamp, double* Geodetic) {
    double t = get_t_second(timestamp);
    get_target_geogetic(id, t, Geodetic);
}

// 获取任意时间戳下的目标地心惯性系坐标（速度未赋值，建议初始化）
inline void get_target_R(const int& id, const std::string& timestamp, double* RV) {
    double t = get_t_second(timestamp);
    get_target_R(id, t, RV);
}

#endif  // ORBIT_PROPAGATION_H
