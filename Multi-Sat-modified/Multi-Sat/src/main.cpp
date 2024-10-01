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

const std::string space = " ";
 
int main() {
	//// ����ʱ���ת��
	//const std::string timestamp = "2035-09-26 12:30:00";
	//double seconds = get_t_second(timestamp);
	//std::cout << "The time t is " << seconds << " sec." << std::endl;

	//// ����ʱ�����ȡĿ���������
	//const int target_id = 1;
	//double geodetic[2];
	//get_target_geogetic(target_id, timestamp, geodetic);
	//std::cout << "Latitude = " << geodetic[0] / D2R << std::endl;
	//std::cout << "Longtitude = " << geodetic[1] / D2R << std::endl;

	//// ���Ծֲ��Ż���к���
	//const double coe0[6] = { 6958.137, 0.000, 2.44, 3.135796327, 3.141592654, 1.3141592654 };
	//double RV0[6]; 
	//int flag;
	//coe2rv(flag, RV0, coe0, mu_km_s);

	//const double next_target_timestamp = get_t_second("2035-09-26 20:00:00");
	//double next_target_geodetic[2], tf, dv[3], RVf[6];
	//get_target_geogetic(target_id, next_target_timestamp, next_target_geodetic);
	//shooting_target2target(10.0, RV0, next_target_timestamp, next_target_geodetic[1], next_target_geodetic[0], tf, dv, RVf, 1);
	//std::cout << "dv[0] = " << std::setprecision(14) << dv[0] << std::endl;
	//std::cout << "dv[1] = " << std::setprecision(14) << dv[1] << std::endl;
	//std::cout << "dv[2] = " << std::setprecision(14) << dv[2] << std::endl;
	//std::cout << "Time stamp = " << std::setprecision(14) << next_target_timestamp << " sec" << std::endl;
	//std::cout << "Real time = " << std::setprecision(14) << tf << " sec" << std::endl;

	MultiTree multi_tree(10000, 4, 50, 1000.0);
	multi_tree.Run();

	


	return 0;
}
