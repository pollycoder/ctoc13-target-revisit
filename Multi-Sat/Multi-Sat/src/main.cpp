/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: main.cpp
* ���ݼ�����GTOC9����������������
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01       2021-05-12    ����      �����ļ�
****************************************************************************/

#include "main.h"

//���ݼ�debris_data
double sats_rv0[5][6];
double city_info[20][2];

/****************************************************************************
* ������   : load_input()
* ��  ��   : ��ȡ���Ǻͳ�����Ϣ��Ԥ������
****************************************************************************/
void load_input()
{
	std::ifstream fin("../input_data/in_SatRV0.txt"); //������Ϣ

	for (int line = 0; line < 5; line++)
	{
		double mass = 1000;
		double rv0[6],coe[6];

		for (int i = 0; i < 6; i++) {
			fin >> rv0[i];// ��λͳһΪkm,s
		}

		//int flag = 0;
		//coe2rv(flag, rv0, coe, mu_km_s);
		memcpy(sats_rv0[line], rv0, 6 * sizeof(double));
	}
	
	fin.close();

	std::ifstream fin2("../input_data/in_PointLonLat.txt"); //������Ϣ

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
	load_input(); //��ȡ�ռ���Ƭ��Ϣ

	MultiTree multi_tree(10000, 4, 50, 300.0);
	multi_tree.Run();

	system("pause");
}
