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

// ����ÿһ�������������֮�ڽ��еĹ۲��ܴ���
int single_sat_score(const double* coe0);
void create_db_single(std::ofstream& fout0);

// ����ÿһ�������������֮�ڽ��еĺ��Ϲ۲��ܴ���
int single_sat_target(const double* coe0, const int& target_id);
void create_db_target(std::ofstream& fout0, const int& target_id);

// ��ȡship�����ݿ�bin��ɸ���۲�������ֵ������
bool filter_bin_file(const std::string& input_filename, const std::string& output_filename, int threshold);

// �����Ǹ������������Ĺ۲�ʱ�̱�
void AccessTableSingleSat(const double* coe0, const std::vector<std::vector<double>>& t_dv_list, std::vector<std::vector<double>>& AccessTable, double& score);

// ���ָ���Ǹ������������Ĺ۲�ʱ�̱�
// ÿ��tuple���������ǵĳ�ʼ��������Ӧ������ʱ�̱�
void AccessTableMultiSat(const std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, std::vector<std::vector<double>>& AccessTable, double& score);


// Ŀ�꺯��
double obj_func(const std::vector<double>& X, std::vector<double>& grad, void* f_data);


// �Ż������Ĳ������֣���ȡ�������Ż���Ϣ
// X���Ż�������4��һ�飬ÿ���Ӧһ��(t, dv[3])
void get_score_info(const std::vector<double>& X, double* f_data, double& score, 
					std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, 
					std::vector<std::vector<double>>& AccessTable);


double obj_func_coelist(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
void get_score_info(const std::vector<double>& X, double* f_data, double& score,
	std::vector<std::vector<double>>& coe_list, std::vector<std::vector<double>>& AccessTable);

double obj_func_4sat(const std::vector<double>& X, std::vector<double>& grad, void* f_data);

void get_score_info_4sat(const std::vector<double>& X, double& score, std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, std::vector<std::vector<double>>& AccessTable);
