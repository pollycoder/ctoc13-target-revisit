#include "single_sat.h"
#include "OrbitFun.h"
#include "visibility_4_targets.h"
#include "Constant.h"
#include <omp.h>
#include <fstream>


int single_sat_score(const double* coe0) {
	double rv0[6];
	const double t0 = 0.0;
	const double tf = 2.0 * 86400.0;
	const double step = 60.0;
	const double wider_angle = 25.0 * D2R;
	int flag;
	coe2rv(flag, rv0, coe0, mu_km_s);

	std::vector<double> rv0_vec(6);
	std::vector<std::vector<double>> rv0_list;
	memcpy(rv0_vec.data(), rv0, 6 * sizeof(double));
	rv0_list.push_back(rv0_vec);
	std::vector<std::vector<double>> result;
	MultiSat_AccessPointObjects(rv0_list, t0, tf, step, 21, result, wider_angle);

	int num = 0;
	for (auto iter = result.begin(); iter != result.end(); iter++) {
		num += iter->size();
	}

	return num;
}


void create_db_single(std::ofstream& fout0) {
	std::vector<std::string> buffer(num_loops);
	auto beforeTime = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Omega_num; i++) {
		double Omega = Omega_start + i * Omega_mesh;
		for (int j = 0; j < inc_num; j++) {
			double inc = inc_start + j * inc_mesh;
			for (int k = 0; k < h_num; k++) {
				double a = Re_km + h_start + k * h_mesh;
				double coe0[6] = { a, 0.0, inc, Omega, 0.0, 0.0 }; // 初始化 coe0 数组
				int index = i * inc_num * h_num + j * h_num + k;
				buffer[index] = std::to_string(coe0[0]) + " " +
					std::to_string(coe0[1]) + " " +
					std::to_string(coe0[2]) + " " +
					std::to_string(coe0[3]) + " " +
					std::to_string(coe0[4]) + " " +
					std::to_string(coe0[5]) + " " +
					std::to_string(single_sat_score(coe0)) + "\n";
			}
		}
	}

	auto afterTime = std::chrono::steady_clock::now();
	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "计算总耗时：" << duration_second << "秒" << std::endl;
	// 16核，平均一条数据0.002s
	// 500 * 200 * 100条数据，预计5.56小时计算完成


	// 写入文件
	//std::ofstream fout0("../output_result/single_sat.bin", std::ios::out | std::ios::binary);
	if (fout0.is_open()) {
		for (const auto& line : buffer) {
			fout0 << line;
		}
		fout0.close();
	}
	else {
		std::cerr << "Unable to open file";
	}

	auto finalTime = std::chrono::steady_clock::now();
	duration_second = std::chrono::duration<double>(finalTime - afterTime).count();
	std::cout << "写入总耗时：" << duration_second << "秒" << std::endl;
	// 预计0.002s写入完成
}


int single_sat_target(const double* coe0, const int& target_id) {
	double rv0[6];
	const double t0 = 0.0;
	const double tf = 2.0 * 86400.0;
	const double step = 60.0;
	const double wider_angle = 25.0 * D2R;
	int flag;
	coe2rv(flag, rv0, coe0, mu_km_s);

	std::vector<double> rv0_vec(6);
	std::vector<std::vector<double>> rv0_list;
	memcpy(rv0_vec.data(), rv0, 6 * sizeof(double));
	rv0_list.push_back(rv0_vec);
	std::vector<std::vector<double>> result;
	AccessPointObjects(rv0, t0, tf, step, target_id, result, wider_angle);


	return result.back().size();
}


void create_db_target(std::ofstream& fout0, const int& target_id) {
	std::vector<std::string> buffer(num_loops);
	auto beforeTime = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Omega_num; i++) {
		double Omega = Omega_start + i * Omega_mesh;
		for (int j = 0; j < inc_num; j++) {
			double inc = inc_start + j * inc_mesh;
			for (int k = 0; k < h_num; k++) {
				double a = Re_km + h_start + k * h_mesh;
				double coe0[6] = { a, 0.0, inc, Omega, 0.0, 0.0 }; // 初始化 coe0 数组
				int index = i * inc_num * h_num + j * h_num + k;
				buffer[index] = std::to_string(coe0[0]) + " " +
					std::to_string(coe0[1]) + " " +
					std::to_string(coe0[2]) + " " +
					std::to_string(coe0[3]) + " " +
					std::to_string(coe0[4]) + " " +
					std::to_string(coe0[5]) + " " +
					std::to_string(single_sat_target(coe0, target_id)) + "\n";
			}
		}
	}

	auto afterTime = std::chrono::steady_clock::now();
	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "计算总耗时：" << duration_second << "秒" << std::endl;
	// 16核，平均一条数据0.002s
	// 500 * 200 * 150条数据，预计8.33小时计算完成


	// 写入文件
	//std::ofstream fout0("../output_result/single_sat.bin", std::ios::out | std::ios::binary);
	if (fout0.is_open()) {
		for (const auto& line : buffer) {
			fout0 << line;
		}
		fout0.close();
	}
	else {
		std::cerr << "Unable to open file";
	}

	auto finalTime = std::chrono::steady_clock::now();
	duration_second = std::chrono::duration<double>(finalTime - afterTime).count();
	std::cout << "写入总耗时：" << duration_second << "秒" << std::endl;
	// 预计0.01s写入完成
}
