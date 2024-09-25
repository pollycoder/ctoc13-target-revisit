#ifndef _RUNGEKUTTAF_H_
#define _RUNGEKUTTAF_H_
#include<iostream>
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"Constant.h"
using namespace std;

//Rugge-Kutta积分
namespace RK{
//	注意：每次调用RKF78进行变步长积分时前都需将RKFiflag设为1，积分成功后其值为2。RKF78返回的整数为变步长积分的步数
//参数说明：
	//Model:被积函数，即yp=f(t,y,para)，其中y为dim维状态量,yp为状态量的一阶导数,para为计算yp所需的一些参数，用数组保存
	//y:dim维状态量的初值，积分成功后，保存最终状态值
	//para:计算yp所需的一些参数，用数组保存
	//tin:积分起始时间，积分成功后，其值变为积分终止时间tout
	//tout:积分终止时间
	//dim:维数
	//RelTol:相对误差
	//AbsTol:绝对误差，与状态量维数一致
	//RKFiflag输入：1全区间开始积分，-1单步长开始积分
	//		输出：2全区间积分成功，-2单步长积分成功，3相对误差太小，4迭代次数超限，5绝对误差太小，6已达到最小允许步长，7计算右函数的次数超限，8不合理的参数，9不可求解
	//RKFWork:一维工作数组14*dim个元素，调用积分器需先分配。比如状态量为6维，则需要开局大小为84的RKFWork。
	//RKFmaxnfe:最大迭代次数，设计一个较大的整数即可，如1000000。
	//fid:保存积分过程的文件，如设计为NULL，即空值，则不保存过程值，可直接从y中得到最终值
int RKF78(void (*Model)(double t, const double* y, double* yp, const double* para), double* y, const double* para,
		  double& tin, double tout, int dim, double RelTol, double* AbsTol, int& RKFiflag, double* RKFWork,
		  int RKFmaxnfe=100000, FILE* fid=NULL);
//wa:2*dim
void RK4(void (*Model)(double t, const double* y, double* yp, const double* para),
			 const double* y, const double* para, double t, double h, double* s, int dim, double* wa);
//wa:3*dim
int RK4(void (*Model)(double t, const double* y, double* yp, const double* para),
			 double* y, const double* para, double& tin, double tout, double h, int dim, double* wa, FILE* fid=NULL);

//RKFWork:一维工作数组7Xdim个元素。RKFiflag,调用RKF45前需设为1，以便进行全区间变步长积分，积分成功后其值为2
int RKF45(void (*Model)(double t, const double* y, double* yp, const double* para), double* y, const double* para,
		  double& tin, double tout, int dim, double RelTol, double* AbsTol, int& RKFiflag, double* RKFWork,
		  int RKFmaxnfe=100000, FILE* fid=NULL);
}

#endif

