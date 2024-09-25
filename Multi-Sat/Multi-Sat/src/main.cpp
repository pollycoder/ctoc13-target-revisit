/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: main.cpp
* 内容简述：GTOC9并行算例的主程序
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-05-12    张众      创建文件
****************************************************************************/

#include "main.h"

//数据集debris_data
double sats_rv0[5][6];
double city_info[20][2];

/****************************************************************************
* 函数名   : load_input()
* 功  能   : 读取卫星和城市信息，预处理部分
****************************************************************************/
void load_input()
{
	std::ifstream fin("../input_data/in_SatRV0.txt"); //卫星信息

	for (int line = 0; line < 5; line++)
	{
		double mass = 1000;
		double rv0[6],coe[6];

		for (int i = 0; i < 6; i++) {
			fin >> rv0[i];// 单位统一为km,s
		}

		//int flag = 0;
		//coe2rv(flag, rv0, coe, mu_km_s);
		memcpy(sats_rv0[line], rv0, 6 * sizeof(double));
	}
	
	fin.close();

	std::ifstream fin2("../input_data/in_PointLonLat.txt"); //城市信息

	for (int line = 0; line < 20; line++)
	{
		double lon, lat;
		fin2 >> lon;
		fin2 >> lat;

		city_info[line][0] = lon * D2R;
		city_info[line][1] = lat * D2R;
	}

	fin2.close();
}
 
int main()
{
	load_input(); //读取空间碎片信息

	MultiTree multi_tree(10000, 4, 50, 300.0);
	multi_tree.Run();

	system("pause");
}
