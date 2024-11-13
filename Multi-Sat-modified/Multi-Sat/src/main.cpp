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
	/*MultiTree multi_tree(10000, 4, 50, 0.3);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "�ܺ�ʱ��" << duration_second << "��" << std::endl;*/
	/***********************************************************************/


	/***********************************************************************/
	// single_imp��ʧ������
	/*double m0 = 1000.0, mf, t0 = 0.0, tf;
	double geo[2];
	get_target_geogetic(2, 0.0, geo);
	double dv[3];
	int NR = 2, branch = 1, flag;
	double coe_00[6] = { sats_coe0[0][0], sats_coe0[0][1], sats_coe0[0][2], sats_coe0[0][3], sats_coe0[0][4], sats_coe0[0][5] };
	double RV0[6];
	coe2rv(flag, RV0, coe_00, mu_km_s);
	single_imp(m0, t0, RV0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, branch);

	std::cout << "flag = " << flag << std::endl;
	std::cout << "tf = " << tf << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;*/


	/***********************************************************************/
	// �����ͬ�������ܹ۲����
	// create_db_single();
	/***********************************************************************/

	/***********************************************************************/
	// �����ͬ�������ض�Ŀ��۲����
	// ��ǰĿ�꣺ship
	// ˳�棺0
	int target_id = 21;
	int branch = 0;
	std::string bra_str;
	if (branch) { bra_str = "R"; }		// Retrograde������
	else { bra_str = "P"; }				// Prograde��˳��
	std::string id = std::to_string(target_id);
	std::string db_filepath = "../output_result/single_sat_" + id + "_" + bra_str + ".bin";
	std::ofstream fout0(db_filepath, std::ios::out | std::ios::binary);
	create_db_target(fout0, target_id);
	// 16�ˣ�ƽ��һ������0.002s
	// 500 * 200 * 150���ݣ�Ԥ��8.33h���
	/***********************************************************************/

	return 0;
}
