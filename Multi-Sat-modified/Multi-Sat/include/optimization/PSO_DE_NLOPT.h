#ifndef _PSO_DE_H_
#define _PSO_DE_H_
#include<vector>
#include "nlopt.hpp"

//PSO�ǲ��а�
//���룺	��С��Ŀ�꺯��
//			1. ObjFun��Ŀ�꺯������������XΪ�Ż�������gradΪ�ݶȣ��������룩��f_dataΪ��ֵ����
//			2. f_dataΪ��ֵ����
//			3. xbestΪ�õ��������Ż�������ÿ����������[0.0, 1.0]֮��
//			4. fbestΪ�õ�������Ŀ����ֵ
//			5. num_variableΪ�Ż������ĸ���
//			6. NpΪPSO�㷨����Ⱥ�ĸ�����Ĭ��������ʱ��Np���㷨���Զ�����Ϊ10*num_variable
//			7. ItMaxΪPSO�㷨������������Ĭ��Ϊ1000
//			8. ItOutΪPSO�㷨��ÿ��ItOut�����һ�����Ž����Ĭ��Ϊ50
void PSO(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000,
	int ItOut = 50, double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5,
	double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.8);

//DE�ǲ��а�
//���룺	��С��Ŀ�꺯��
//			1. ObjFun��Ŀ�꺯������������XΪ�Ż�������gradΪ�ݶȣ��������룩��f_dataΪ��ֵ����
//			2. f_dataΪ��ֵ����
//			3. xbestΪ�õ��������Ż�������ÿ����������[0.0, 1.0]֮��
//			4. fbestΪ�õ�������Ŀ����ֵ
//			5. num_variableΪ�Ż������ĸ���
//			6. NpΪDE�㷨����Ⱥ�ĸ�����Ĭ��������ʱ��Np���㷨���Զ�����Ϊ10*num_variable
//			7. ItMaxΪDE�㷨������������Ĭ��Ϊ1000
//			8. ItOutΪDE�㷨��ÿ��ItOut�����һ�����Ž����Ĭ��Ϊ50
//			9. CRΪDE�㷨�Ľ����ʣ�Ĭ��Ϊ-1ʱ��CR���㷨���Զ�����Ϊmin((18.0 + num_variable) / 40.0, 0.9)
void DE(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000, int ItOut = 50, double CR = -1.0);

//PSO���а�
//���룺	��С��Ŀ�꺯��
//			1. ObjFun��Ŀ�꺯������������XΪ�Ż�������gradΪ�ݶȣ��������룩��f_dataΪ��ֵ����
//			2. f_dataΪ��ֵ����
//			3. xbestΪ�õ��������Ż�������ÿ����������[0.0, 1.0]֮��
//			4. fbestΪ�õ�������Ŀ����ֵ
//			5. num_variableΪ�Ż������ĸ���
//			6. NpΪPSO�㷨����Ⱥ�ĸ�����Ĭ��������ʱ��Np���㷨���Զ�����Ϊ10*num_variable
//			7. ItMaxΪPSO�㷨������������Ĭ��Ϊ1000
//			8. ItOutΪPSO�㷨��ÿ��ItOut�����һ�����Ž����Ĭ��Ϊ50
void PSO_parallel(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000,
	int ItOut = 50, double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5,
	double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.8);

//DE���а�
//���룺	��С��Ŀ�꺯��
//			1. ObjFun��Ŀ�꺯������������XΪ�Ż�������gradΪ�ݶȣ��������룩��f_dataΪ��ֵ����
//			2. f_dataΪ��ֵ����
//			3. xbestΪ�õ��������Ż�������ÿ����������[0.0, 1.0]֮��
//			4. fbestΪ�õ�������Ŀ����ֵ
//			5. num_variableΪ�Ż������ĸ���
//			6. NpΪDE�㷨����Ⱥ�ĸ�����Ĭ��������ʱ��Np���㷨���Զ�����Ϊ10*num_variable
//			7. ItMaxΪDE�㷨������������Ĭ��Ϊ1000
//			8. ItOutΪDE�㷨��ÿ��ItOut�����һ�����Ž����Ĭ��Ϊ50
//			9. CRΪDE�㷨�Ľ����ʣ�Ĭ��Ϊ-1ʱ��CR���㷨���Զ�����Ϊmin((18.0 + num_variable) / 40.0, 0.9)
void DE_parallel(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000, int ItOut = 50, double CR = -1.0);

void nlopt_main(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
	void* f_data, std::vector<double>& X, double& f, int num_variable, int Np = 0, int maxstep = 1000);
#endif
