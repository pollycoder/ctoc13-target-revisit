#pragma once
#include "Constant.h"
#include <vector>
#include <string>
#include <tuple>

// Height: 800km - 1000km
// RAAN: 0 - 360 deg
// INC: 40 - 60 deg
const double h_start = 800.0; const double h_end = 1000.0;
const double Omega_start = 0.0; const double Omega_end = D2PI;
const double inc_start = 40.0 * D2R; const double inc_end = 60.0 * D2R;

// Mesh: 500 * 200 * 150
// 8.33h for each database
const int Omega_num = 500;
const int inc_num = 200;
const int h_num = 150;
const int num_loops = Omega_num * inc_num * h_num;
const double h_mesh = (h_end - h_start) / h_num;
const double Omega_mesh = (Omega_end - Omega_start) / Omega_num;
const double inc_mesh = (inc_end - inc_start) / inc_num;

// 计算每一个轨道根数两天之内进行的观测总次数
int single_sat_score(const double* coe0);
void create_db_single(std::ofstream& fout0);

// 计算每一个轨道根数两天之内进行的海上观测总次数
int single_sat_target(const double* coe0, const int& target_id);
void create_db_target(std::ofstream& fout0, const int& target_id);

// 读取ship的数据库bin，筛掉观测少于阈值的数据
bool filter_bin_file(const std::string& input_filename, const std::string& output_filename, int threshold);

// 单颗星给定机动次数的观测时刻表
void AccessTableSingleSat(const double* coe0, const std::vector<std::vector<double>>& t_dv_list, std::vector<std::vector<double>>& AccessTable, double& score);

// 多颗指定星给定机动次数的观测时刻表
// 每个tuple：单颗卫星的初始根数，对应的脉冲时刻表
void AccessTableMultiSat(const std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, std::vector<std::vector<double>>& AccessTable, double& score);


// 目标函数
double obj_func(const std::vector<double>& X, std::vector<double>& grad, void* f_data);


// 优化函数的操作部分：获取分数和优化信息
// X：优化变量，4个一组，每组对应一个(t, dv[3])
void get_score_info(const std::vector<double>& X, double* f_data, double& score, 
					std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, 
					std::vector<std::vector<double>>& AccessTable);


double obj_func_coelist(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
void get_score_info(const std::vector<double>& X, double* f_data, double& score,
	std::vector<std::vector<double>>& coe_list, std::vector<std::vector<double>>& AccessTable);

double obj_func_4sat(const std::vector<double>& X, std::vector<double>& grad, void* f_data);

void get_score_info_4sat(const std::vector<double>& X, double& score, std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, std::vector<std::vector<double>>& AccessTable);
