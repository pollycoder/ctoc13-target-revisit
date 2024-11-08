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

 
int main() {
	/***********************************************************************/
	// ����
	/*MultiTree multi_tree(10, 4, 50, 0.3);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "�ܺ�ʱ��" << duration_second << "��" << std::endl;*/


	/***********************************************************************/
	// single_imp��ʧ������
	/*double m0 = 1000.0, mf, t0 = 0.0, tf;
	double geo[2];
	get_target_geogetic(2, 0.0, geo);
	double dv[3];
	int NR = 1, branch = 0, flag;
	double coe_00[6] = { sats_coe0[1][0], sats_coe0[1][1], sats_coe0[1][2], sats_coe0[1][3], sats_coe0[1][4], sats_coe0[1][5] };
	double RV0[6];
	coe2rv(flag, RV0, coe_00, mu_km_s);
	single_imp(m0, t0, RV0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, branch);

	std::cout << "flag = " << flag << std::endl;
	std::cout << "tf = " << tf << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;*/


	/***********************************************************************/
	// �����ͬ�����Ĺ۲����
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
	std::ofstream fout0("../output_result/single_sat.bin", std::ios::out | std::ios::binary);
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

	return 0;
}
