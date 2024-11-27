#include "single_sat.h"
#include "OrbitFun.h"
#include "visibility_4_targets.h"
#include "Constant.h"
#include <omp.h>
#include <fstream>
#include "max_reseetime.h"
#include<numeric>
#include <algorithm>


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


// 读取ship的数据库bin，筛掉观测少于阈值的数据
bool filter_bin_file(const std::string& input_filename, const std::string& output_filename, int threshold) {
	// 打开输入二进制文件
	std::ifstream infile(input_filename);
	if (!infile) {
		std::cerr << "无法打开输入文件: " << input_filename << std::endl;
		return false;
	}

	// 存储符合条件的数据
	std::vector<std::vector<double>> filtered_lines;

	while (true) {
		// 每行数据
		std::vector<double> line(7);
		for (int i = 0; i < 7; ++i) {
			if (!(infile >> line[i])) {
				// 到文件末尾，退出循环
				break;
			}
		}

		// 如果读取到的数据不足7个，则认为是文件结束
		if (line.size() < 7) {
			break;
		}

		double seventh_element = line[6];
		if (seventh_element >= threshold) {
			filtered_lines.push_back(line);
		}
	}

	// 关闭输入文件
	infile.close();

	// 打开输出二进制文件
	std::ofstream outfile(output_filename);
	if (!outfile) {
		std::cerr << "无法打开输出文件: " << output_filename << std::endl;
		return false;
	}

	for (const auto& line : filtered_lines) {
		for (int i = 0; i < 7; ++i) {
			outfile << line[i] << ' ';
		}
		outfile << '\n';
	}

	// 关闭输出文件
	outfile.close();

	return true;
}

void AccessTableSingleSat(const double* coe0, const std::vector<std::vector<double>>& t_dv_list, std::vector<std::vector<double>>& AccessTable, double& score) {
	int flag;
	double rv0[6], rv0_temp[6], rv_imp[6];															// 初始位置速度，脉冲点位置速度
	std::vector<std::vector<double>> timetable_temp;

	AccessTable.clear();
	AccessTable.resize(TargetNum);


	coe2rv(flag, rv0, coe0, mu_km_s);						
	memcpy(rv0_temp, rv0, 6 * sizeof(double));
	if (t_dv_list.empty()) {
		AccessPointObjects(rv0, 0.0, 172800.0, 60.0, TargetNum, AccessTable, 20.0 * D2R);
		return;
	}

	double t0 = 0, tf = t_dv_list[0][0];
	//std::cout << coe0[0] << " " << coe0[1] << " " << coe0[2] << " " << coe0[3] << " " << coe0[4] << " " << coe0[5] << std::endl;

	// 求解流程：rv0 -> time table -> rvf -> dv
	for (int i = 0; i < t_dv_list.size() + 1; i++) {
		timetable_temp.clear();
		AccessPointObjects(rv0_temp, t0, tf, 60.0, 21, timetable_temp, 20.0 * D2R);
		/*for(const auto& t:timetable_temp[0]) std::cout << t << " ";
		std::cout << std::endl;*/
		for (int j = 0; j < TargetNum; j++) {
			AccessTable[j].insert(AccessTable[j].end(), timetable_temp[j].begin(), timetable_temp[j].end());
		}

		//std::cout << "t0 = " << t0 << " tf = " << tf << std::endl;
		propagate_j2(rv0_temp, rv_imp, t0, tf);
		t0 = tf;
		if (i == t_dv_list.size()) {
			return;
		}
		else {
			if (i == t_dv_list.size() - 1) tf = 172800.0;
			else tf = t_dv_list[i + 1][0];
			for (int j = 0; j < 3; j++) rv_imp[j + 3] += t_dv_list[i][j + 1];
			double coe[6];
			rv2coe(flag, coe, rv_imp, mu_km_s);
			const double apo = coe[0] * (1 + coe[1]) - Re_km;
			const double peri = coe[0] * (1 - coe[1]) - Re_km;
			score += (std::max(apo, hmax) - hmax) * (std::max(apo, hmax) - hmax);
			score += (std::min(peri, hmin) - hmin) * (std::min(peri, hmin) - hmin);
		}

		memcpy(rv0_temp, rv_imp, 6 * sizeof(double));
	}
}


void AccessTableMultiSat(const std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, std::vector<std::vector<double>>& AccessTable, double& score) {
	AccessTable.clear();
	AccessTable.resize(TargetNum);
	for (const auto& sat : sat_info_list) {
		std::vector<std::vector<double>> accesstable_temp;
		const std::vector<double> coe0_vec = std::get<0>(sat);
		const std::vector<std::vector<double>> t_dv_list = std::get<1>(sat);
		const double coe0[6] = { coe0_vec[0], coe0_vec[1], coe0_vec[2], coe0_vec[3], coe0_vec[4], coe0_vec[5] };
		AccessTableSingleSat(coe0, t_dv_list, accesstable_temp, score);
		for (int i = 0; i < TargetNum; i++) {
			AccessTable[i].insert(AccessTable[i].end(), accesstable_temp[i].begin(), accesstable_temp[i].end());
		}
	}

	for (auto& row : AccessTable) {
		std::sort(row.begin(), row.end());
	}
}

double obj_func(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	//auto beforeTime = std::chrono::steady_clock::now();
	double* para = static_cast<double*>(f_data);
	double score;
	std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>> sat_info_list;
	std::vector<std::vector<double>> AccessTable;

	get_score_info(X, para, score, sat_info_list, AccessTable);
	// afterTime = std::chrono::steady_clock::now();
	//double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	//std::cout << duration_second << std::endl;
	return score;
}


void get_one_sat(int& i, const std::vector<double>& X, int& imp, int& imp_idx, double& score, std::vector<std::vector<double>>& t_dv, std::vector<double>& coe0, std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list)
{
	int fixed_idx = 0;
	// 如果是脉冲星,按照次数加入脉冲信息，不是则信息表为空
	if (std::find(imp_sat.begin(), imp_sat.end(), i) != imp_sat.end()) {
		double t_temp = 0.0;
		for (int j = 0; j < imp_num[imp_idx]; j++) {
			t_temp += X[imp * 4] * 172800.0;
			score += (std::max(t_temp, 172800.0) - 172800.0) * (std::max(t_temp, 172800.0) - 172800.0);
			double dv[3] = {
				(X[imp * 4 + 1] - 0.5) * imp_max,
				(X[imp * 4 + 2] - 0.5) * imp_max,
				(X[imp * 4 + 3] - 0.5) * imp_max
			};
			std::vector<double> t_dv_one = { t_temp, dv[0], dv[1], dv[2] };
			t_dv.push_back(t_dv_one);
			imp++;
		}
		imp_idx++;
	} else {
		for (const auto& idx : fixed_sat) {
			if (idx == i) {
				t_dv = fixed_imp[fixed_idx];
				break;
			}
			else {
				fixed_idx++;
			}
		}
	}

	std::tuple<std::vector<double>, std::vector<std::vector<double>>> sat = std::make_tuple(coe0, t_dv);
	sat_info_list.push_back(sat);
}

void get_score_info(const std::vector<double>& X, double* f_data, double& score,
	std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list,
	std::vector<std::vector<double>>& AccessTable) {
	score = 0.0;

	int imp = 0;											// 记录已加脉冲的次数
	int imp_idx = 0;										// 记录已加脉冲星的个数
	for (int i = 0; i < TreeNum; i++) {
		std::vector<double> coe0;
		for (int j = 0; j < 6; j++) coe0.push_back(sats_coe0[i][j]);
		std::vector<std::vector<double>> t_dv;

		get_one_sat(i, X, imp, imp_idx, score, t_dv, coe0, sat_info_list);
	}

	AccessTableMultiSat(sat_info_list, AccessTable, score);
	
	std::vector<double> max_revisit_gap;
	max_reseetime(AccessTable, max_revisit_gap);
	int idx = 0;
	for (const auto& gap : max_revisit_gap) {
		idx++;
		if (idx != TargetNum) 
			score -= 4.0 * 6.0 / std::max(6.0, gap);
		else
			score -= 20.0 * 3.0 / std::max(3.0, gap);
		
	}

}