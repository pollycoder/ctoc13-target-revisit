#ifndef _PSO_DE_H_
#define _PSO_DE_H_
#include<vector>
#include "nlopt.hpp"

//PSO非并行版
//输入：	最小化目标函数
//			1. ObjFun，目标函数。其输入量X为优化变量，grad为梯度（不必输入），f_data为常值参数
//			2. f_data为常值阐述
//			3. xbest为得到的最优优化变量，每个分量均在[0.0, 1.0]之间
//			4. fbest为得到的最优目标数值
//			5. num_variable为优化变量的个数
//			6. Np为PSO算法中种群的个数，默认输入零时，Np在算法中自动设置为10*num_variable
//			7. ItMax为PSO算法最大迭代次数，默认为1000
//			8. ItOut为PSO算法中每隔ItOut代输出一次最优结果，默认为50
void PSO(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000,
	int ItOut = 50, double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5,
	double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.8);

//DE非并行版
//输入：	最小化目标函数
//			1. ObjFun，目标函数。其输入量X为优化变量，grad为梯度（不必输入），f_data为常值参数
//			2. f_data为常值阐述
//			3. xbest为得到的最优优化变量，每个分量均在[0.0, 1.0]之间
//			4. fbest为得到的最优目标数值
//			5. num_variable为优化变量的个数
//			6. Np为DE算法中种群的个数，默认输入零时，Np在算法中自动设置为10*num_variable
//			7. ItMax为DE算法最大迭代次数，默认为1000
//			8. ItOut为DE算法中每隔ItOut代输出一次最优结果，默认为50
//			9. CR为DE算法的交叉率，默认为-1时，CR在算法中自动设置为min((18.0 + num_variable) / 40.0, 0.9)
void DE(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000, int ItOut = 50, double CR = -1.0);

//PSO并行版
//输入：	最小化目标函数
//			1. ObjFun，目标函数。其输入量X为优化变量，grad为梯度（不必输入），f_data为常值参数
//			2. f_data为常值阐述
//			3. xbest为得到的最优优化变量，每个分量均在[0.0, 1.0]之间
//			4. fbest为得到的最优目标数值
//			5. num_variable为优化变量的个数
//			6. Np为PSO算法中种群的个数，默认输入零时，Np在算法中自动设置为10*num_variable
//			7. ItMax为PSO算法最大迭代次数，默认为1000
//			8. ItOut为PSO算法中每隔ItOut代输出一次最优结果，默认为50
void PSO_parallel(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000,
	int ItOut = 50, double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5,
	double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.8);

//DE并行版
//输入：	最小化目标函数
//			1. ObjFun，目标函数。其输入量X为优化变量，grad为梯度（不必输入），f_data为常值参数
//			2. f_data为常值阐述
//			3. xbest为得到的最优优化变量，每个分量均在[0.0, 1.0]之间
//			4. fbest为得到的最优目标数值
//			5. num_variable为优化变量的个数
//			6. Np为DE算法中种群的个数，默认输入零时，Np在算法中自动设置为10*num_variable
//			7. ItMax为DE算法最大迭代次数，默认为1000
//			8. ItOut为DE算法中每隔ItOut代输出一次最优结果，默认为50
//			9. CR为DE算法的交叉率，默认为-1时，CR在算法中自动设置为min((18.0 + num_variable) / 40.0, 0.9)
void DE_parallel(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
	std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000, int ItOut = 50, double CR = -1.0);

void nlopt_main(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
	void* f_data, std::vector<double>& X, double& f, int num_variable, int Np = 0, int maxstep = 1000);
#endif
