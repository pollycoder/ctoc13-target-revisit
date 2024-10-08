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
****************************************************************************/
#include <stdlib.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"

const std::string space = " ";
 
int main() {
	MultiTree multi_tree(100, 3, 10, 1000.0);
	multi_tree.Run();

	/*std::vector<double> X = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double f = 1000.0;
	DE_parallel(obj_func_initial_coe, nullptr, X, f, 6, 1000, 10000, 100);

	double coe0[6];
	get_initial_coe(X, coe0, f);
	std::cout << coe0[0] << space << coe0[1] << space << coe0[2] << space << coe0[3] << space << coe0[4] << space << coe0[5] << std::endl;*/

	test_lambert();


	return 0;
}
