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
				double coe0[6] = { a, 0.0, inc, Omega, 0.0, 0.0 }; // ��ʼ�� coe0 ����
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
	std::cout << "�����ܺ�ʱ��" << duration_second << "��" << std::endl;
	// 16�ˣ�ƽ��һ������0.002s
	// 500 * 200 * 100�����ݣ�Ԥ��5.56Сʱ�������


	// д���ļ�
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
	std::cout << "д���ܺ�ʱ��" << duration_second << "��" << std::endl;
	// Ԥ��0.002sд�����
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
				double coe0[6] = { a, 0.0, inc, Omega, 0.0, 0.0 }; // ��ʼ�� coe0 ����
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
	std::cout << "�����ܺ�ʱ��" << duration_second << "��" << std::endl;
	// 16�ˣ�ƽ��һ������0.002s
	// 500 * 200 * 150�����ݣ�Ԥ��8.33Сʱ�������


	// д���ļ�
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
	std::cout << "д���ܺ�ʱ��" << duration_second << "��" << std::endl;
	// Ԥ��0.01sд�����
}


// ��ȡship�����ݿ�bin��ɸ���۲�������ֵ������
bool filter_bin_file(const std::string& input_filename, const std::string& output_filename, int threshold) {
	// ������������ļ�
	std::ifstream infile(input_filename);
	if (!infile) {
		std::cerr << "�޷��������ļ�: " << input_filename << std::endl;
		return false;
	}

	// �洢��������������
	std::vector<std::vector<double>> filtered_lines;

	while (true) {
		// ÿ������
		std::vector<double> line(7);
		for (int i = 0; i < 7; ++i) {
			if (!(infile >> line[i])) {
				// ���ļ�ĩβ���˳�ѭ��
				break;
			}
		}

		// �����ȡ�������ݲ���7��������Ϊ���ļ�����
		if (line.size() < 7) {
			break;
		}

		double seventh_element = line[6];
		if (seventh_element >= threshold) {
			filtered_lines.push_back(line);
		}
	}

	// �ر������ļ�
	infile.close();

	// ������������ļ�
	std::ofstream outfile(output_filename);
	if (!outfile) {
		std::cerr << "�޷�������ļ�: " << output_filename << std::endl;
		return false;
	}

	for (const auto& line : filtered_lines) {
		for (int i = 0; i < 7; ++i) {
			outfile << line[i] << ' ';
		}
		outfile << '\n';
	}

	// �ر�����ļ�
	outfile.close();

	return true;
}


// Tool functions (static)
static void pert_single_sat(const std::vector<double>& X, const double* para, double& t, double* dv, double* coe0) {
	for (int i = 0; i < 6; i++) coe0[i] = para[i];
	const double dt = 43200.0;										// ʱ������Ŷ���ǰ��6Сʱ
	const double d_imp = 1.0;										// �����������Ŷ�������0.5km/s
	t = para[6] + (X[0] - 0.5) * dt;
	dv[0] = (X[1] - 0.5) * d_imp;
	dv[1] = (X[2] - 0.5) * d_imp;
	dv[2] = (X[3] - 0.5) * d_imp;
}

static void get_imp_trajectory(int& flag, double* rv0, const double* coe0, double* rv_imp, double& t_imp, double* dv, double& peri, double& apo)
{
	coe2rv(flag, rv0, coe0, mu_km_s);
	propagate_j2(rv0, rv_imp, 0.0, t_imp, 1.0e-5, 1.0e-7);
	for (int i = 0; i < 3; i++) rv_imp[i + 3] += dv[i];
	double coe_imp[6], a, e;
	rv2coe(flag, coe_imp, rv_imp, mu_km_s);
	a = coe_imp[0];
	e = coe_imp[1];
	peri = a * (1 - e) - Re_km;
	apo = a * (1 + e) - Re_km;
}

double obj_single_sat(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	const double* para = static_cast<double*>(f_data);
	double t_impulse, dv[3], mean;
	std::vector<double> max_revisit;
	get_revisit(X, para, max_revisit, t_impulse, dv, mean);

	return mean;
}

void get_revisit(const std::vector<double>& X, const double* para, std::vector<double>& max_revisit, double& t_imp, double* dv, double& score) {
	double coe0[6], rv0[6], rv_imp[6], peri, apo;
	pert_single_sat(X, para, t_imp, dv, coe0);
	int flag;
	get_imp_trajectory(flag, rv0, coe0, rv_imp, t_imp, dv, peri, apo);
	if (peri < 201.0 || apo > 999.0) {
		score =  penalty;
		return;
	}

	std::vector <std::vector<double>> visible_list;
	std::vector<std::vector<double>> append_list;
	AccessPointObjects(rv0, 0.0, t_imp, 60.0, 21, visible_list);
	AccessPointObjects(rv_imp, t_imp, 2.0 * 86400.0, 60.0, 21, append_list);
	for (int i = 0; i < 21; i++) {
		visible_list[i].insert(visible_list[i].end(), append_list[i].begin(), append_list[i].end());
	}

	max_reseetime(visible_list, max_revisit);

	double sum = std::accumulate(max_revisit.begin(), max_revisit.end(), 0.0);
	score = sum / 21.0;
	
}


// �Ż�����ǵ��λ���
// ���������7��һѭ������
//		f_data[0-5]���������
//		f_data[6]������ʱ��
// �Ż�������4��һѭ������
//		t������ʱ�䣬��ֵ�ݶ�Ϊ1�죬�Ŷ���Χ6Сʱ
//		dv[3]������3����������ֵΪ0���Ŷ���Χ0.5km/s
// Ŀ�꣺
//		t_gap_ave��ƽ��ÿ��Ŀ����ط�ʱ��
// Լ����
//		h������߶ȣ�200km-1000km
//		dv�����岻����1km/s���Ŷ���Χע�����ᷢ�������⣬���ԣ�
double obj_multi_sat(const std::vector<double>& X, std::vector<double>& grad, void* f_data) {
	const double* para = static_cast<double*>(f_data);
	double f;
	std::vector<std::vector<double>> dv_list;
	std::vector<double> t_imp_list, max_revisit;
	get_revisit(X, para, max_revisit, t_imp_list, dv_list, f);
	return f;
}

void get_revisit(const std::vector<double>& X, const double* para, std::vector<double>& max_revisit, std::vector<double>& t_imp, std::vector<std::vector<double>>& dv, double& score) {
	const int paranum_group = 7;
	const int varnum_group = 4;
	
	std::vector<std::vector<double>> visible_list;
	visible_list.resize(21);
	int flag;

	for (int i = 0; i < TreeNum; i++) {
		std::vector<double> X_current;
		for (int j = 0; j < varnum_group; j++) X_current.push_back(X[i * varnum_group + j]);

		double coe0_current[6], rv0_current[6], t_imp_current,
			   dv_current[3], rv_imp_current[6], peri, apo;
		double para_current[paranum_group];
		for (int j = 0; j < paranum_group; j++) para_current[j] = para[i * paranum_group + j];

		pert_single_sat(X_current, para_current, t_imp_current, dv_current, coe0_current);

		get_imp_trajectory(flag, rv0_current, coe0_current, rv_imp_current, t_imp_current, dv_current, peri, apo);
		if (peri < 201.0 || apo > 999.0) {
			score = penalty;
			return;
		}

		const std::vector<double> dv_current_vec = { dv_current[0], dv_current[1], dv_current[2] };
		std::vector <std::vector<double>> visible_list_current;
		std::vector<std::vector<double>> append_list_current;
		AccessPointObjects(rv0_current, 0.0, t_imp_current, 60.0, 21, visible_list_current);
		AccessPointObjects(rv_imp_current, t_imp_current, 2.0 * 86400.0, 60.0, 21, append_list_current);
		for (int i = 0; i < 21; i++) {
			visible_list[i].insert(visible_list[i].end(), visible_list_current[i].begin(), visible_list_current[i].end());
			visible_list[i].insert(visible_list[i].end(), append_list_current[i].begin(), append_list_current[i].end());
		}

		t_imp.push_back(t_imp_current);
		dv.push_back(dv_current_vec);
	}

	for (int i = 0; i < 21; i++) {
		std::sort(visible_list[i].begin(), visible_list[i].end());
	}

	max_reseetime(visible_list, max_revisit);
	score = std::accumulate(max_revisit.begin(), max_revisit.end(), 0.0);
	score /= 21.0;
}