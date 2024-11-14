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


// ����������������⣬���ã�
void multitree_search() {
	MultiTree multi_tree(1000, 4, 50, 0.3);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "�ܺ�ʱ��" << duration_second << "��" << std::endl;
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
	get_target_geogetic(0, 0.0, geo);
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

 
int main() {
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

	return 0;
}
