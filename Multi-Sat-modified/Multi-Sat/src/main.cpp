/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: main.cpp
* ���ݼ�����GTOC9���������������򣬺��޸�����CTOC13����
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01       2021-05-12    ����      �����ļ�
* 02	   2024-09-28	 ��ܲ	   ���Թ۲��к���
* 03	   2024-09-30	 ��		   ����������
* 04	   2024-10-15	 ��		   ��ʼ��ֵ̽��
****************************************************************************/
#include <stdlib.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"
#include "J2propagation.h"
#include "single_sat.h"
#include <iomanip>



const std::string space = " ";


// �����������������
void multitree_search() {
	MultiTree multi_tree(1000, 4, 50, dv_max);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "�ܺ�ʱ��" << duration_second << "��" << std::endl;
	//0.01s��չһ���ڵ�
	//W=10000����չһ��130s
	//����176��Լ22880s��Ԥ��7h����
	//������չԼ300-400�㣬һ���ڱ������
}


/***********************************************************************/
// single_imp��ʧ������
// Ŀ�꣺12
// Ԥ�ƹ۲�ʱ�䣺1620s
// ���⣺tG���ƴ���84000s������Ȧ�����ز�ƥ��ʱ�����޽�
// Ŀ�꣺13
// Ԥ�ƹ۲�ʱ�䣺1800s
// ���⣺��������ȷ����Ȧ����ƥ�䣨NR=3��ʵ��1Ȧ���ܹ۲⣩��������ܴ�dv = 0.11047 -3.51521 -2.95107��
// Ŀ�꣺14
// Ԥ�ƹ۲�ʱ�䣺8400s������ǳ���ȷ��
// ���⣺��������ȷ����Ȧ����ƥ�䣨NR=3��ʵ��2Ȧ���ܹ۲⣩�������岻�ӽ�0��dv = 0.00431052 -0.137163 -0.11515��
// ��13000s�ڶ��ι۲�Ŀ��12��ʼ�����������ȷ
/***********************************************************************/
void fail_example_singleimp() {
	double m0 = 1000.0, mf, t0 = 0.0, tf;
	double geo[2];
	get_target_geogetic(11, 0.0, geo);
	double dv[3];
	int NR = 4, branch = 0, flag;
	double coe_00[6] = { sats_coe0[0][0], sats_coe0[0][1], sats_coe0[0][2], sats_coe0[0][3], sats_coe0[0][4], sats_coe0[0][5] };
	double RV0[6];
	coe2rv(flag, RV0, coe_00, mu_km_s);
	single_imp(m0, t0, RV0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, branch);

	std::cout << "flag = " << flag << std::endl;
	std::cout << "tf = " << tf << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
}


// ��ȡ���ݿ⣬ɸѡ�������µ����ݿ�
void read_db() {
	const std::string input_filename1 = "../output_result/single_sat_21_R.bin";
	const std::string output_filename1 = "../output_result/single_better_21_R.bin";
	int threshold = 10;

	auto beforeTime = std::chrono::steady_clock::now();
	if (filter_bin_file(input_filename1, output_filename1, threshold)) {
		auto afterTime = std::chrono::steady_clock::now();
		double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
		std::cout << "1�����ݿ⴦����ɣ��ܺ�ʱ��" << duration_second << "��" << std::endl;
	}
	else {
		std::cerr << "����ʧ��" << std::endl;
	}
	// 0.001s��ȡһ�����ݣ�Ԥ��4.16h����һ�����ݿ�

	const std::string input_filename2 = "../output_result/single_sat_21_P.bin";
	const std::string output_filename2 = "../output_result/single_better_21_P.bin";

	beforeTime = std::chrono::steady_clock::now();
	if (filter_bin_file(input_filename2, output_filename2, threshold)) {
		auto afterTime = std::chrono::steady_clock::now();
		double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
		std::cout << "2�����ݿ⴦����ɣ��ܺ�ʱ��" << duration_second << "��" << std::endl;
	}
	else {
		std::cerr << "����ʧ��" << std::endl;
	}
	// 0.001s��ȡһ�����ݣ�Ԥ��4.16h����һ�����ݿ�

	// 16�ˣ��ܼ�8.33h���

}

void single_sat_opt() {
	double para[7] = { sats_coe0[1][0], sats_coe0[1][1], sats_coe0[1][2], sats_coe0[1][3], sats_coe0[1][4], sats_coe0[1][5], 86400.0 };
	std::vector<double> X = { 0.5, 0.5, 0.5, 0.5 };
	double f;

	// DE global
	auto beforeTime = std::chrono::steady_clock::now();
	DE_parallel(obj_single_sat, para, X, f, X.size(), 10, 5000, 100);
	auto afterTime = std::chrono::steady_clock::now();
	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "DE�ܺ�ʱ��" << duration_second << "��" << std::endl;

	std::vector<double> max_revisit;
	double t_imp, dv[3];
	get_revisit(X, para, max_revisit, t_imp, dv, f);

	std::ofstream fout0("../output_result/result.txt");
	std::cout << "t = " << t_imp << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;

	fout0 << sats_coe0[0][0] << " " << sats_coe0[0][1] << " " << sats_coe0[0][2] << " " << sats_coe0[0][3] << " " << sats_coe0[0][4] << " " << sats_coe0[0][5] << std::endl;
	fout0 << t_imp << " " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
	int index = 0;
	for (auto iter = max_revisit.begin(); iter != max_revisit.end(); iter++) {
		index++;
		fout0 << "Target" << index << ": " << *iter << std::endl;
		std::cout << *iter << std::endl;
	}

	// ��Ⱥ��50
	// ����������10000
	// 0.025sĿ�꺯������һ��
	// 16�ˣ�Ԥ��3.33h���
}

void multi_sat_opt() {
	int para[impNum] = {4, 5, 6, 7};
	std::vector<double> X;
	std::vector<int> imp_idx;
	for (int i = 0; i < impNum; i++) {
		for (int j = 0; j < 4; j++) X.push_back(0.5);
		imp_idx.push_back(para[i]);
	}

	double f;

	// DE global
	/*auto beforeTime = std::chrono::steady_clock::now();
	DE_parallel(obj_multi_sat_certain, para, X, f, X.size(), 40, 10000, 50);
	nlopt_main(obj_multi_sat_certain, para, X, f, X.size(), 5000);
	auto afterTime = std::chrono::steady_clock::now();
	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "�ܺ�ʱ��" << duration_second << "��" << std::endl;*/

	std::vector<double> max_revisit;
	std::vector<double> t_imp;
	std::vector<std::vector<double>> dv;
	
	//DE_parallel(obj_func_imp_coe, para, X, f, X.size(), 5 * X.size(), 10000, 100);
	get_revisit_certain(X, para, max_revisit, t_imp, dv, f);
	std::cout << "f = " << f << std::endl;

	std::ofstream fout0("../output_result/result.txt");
	int idx = 0;
	for (int i = 0; i < TreeNum; i++) {
		fout0 << sats_coe0[i][0] << " " << sats_coe0[i][1] << " " << sats_coe0[i][2] << " " << sats_coe0[i][3] << " " << sats_coe0[i][4] << " " << sats_coe0[i][5] << std::endl;
		if (std::find(imp_idx.begin(), imp_idx.end(), i) != imp_idx.end()) {
			fout0 << t_imp[idx] << " " << dv[idx][0] << " " << dv[idx][1] << " " << dv[idx][2] << std::endl;
			idx++;
		}
	}

	int index = 0;
	for (auto iter = max_revisit.begin(); iter != max_revisit.end(); iter++) {
		index++;
		fout0 << "Target" << index << ": " << *iter << std::endl;
		std::cout << "Target" << index << ": " << *iter << std::endl;
	}

	// ��Ⱥ��160
	// ����������10000
	// 0.016sĿ�꺯������һ��
	// 16�ˣ�Ԥ��8h���
}

void max_revisit_verify() {
	std::vector<std::vector<double>> rv0_list;
	std::vector<std::vector<double>> result;
	for (int i = 0; i < TreeNum; i++) {
		std::vector<double> rv0;
		double rv0temp[6];
		int flag;
		coe2rv(flag, rv0temp, sats_coe0[i], mu_km_s);
		for (int j = 0; j < 6; j++) rv0.push_back(rv0temp[j]);
		rv0_list.push_back(rv0);
	}
	MultiSat_AccessPointObjects(rv0_list, 0.0, 172800.0, 60.0, 21, result);
	std::vector<double> max_revisit;
	max_reseetime(result, max_revisit);
	for (auto row = result.begin(); row != result.end(); row++) {
		for (auto it = row->begin(); it != row->end(); it++) {
			std::cout << *it << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
	std::cout << "Max revisit gap:" << std::endl;
	for (auto it = max_revisit.begin(); it != max_revisit.end(); it++) {
		std::cout << *it << std::endl;
	}
}
 
int main() {

	// ������16��
	// ����
	//	1. ��������
	// ��ʼʱ�䣺15:00
	// ��һ������ʱ�䣺21:00
	multitree_search();
	//max_revisit_verify();
	//multi_sat_opt();
	return 0;
}
