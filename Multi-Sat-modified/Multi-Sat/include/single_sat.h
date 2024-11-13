#pragma once
#include "Constant.h"
#include <vector>
#include <string>

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
