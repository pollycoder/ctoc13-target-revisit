/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: model.h
* 内容简述：用于模型信息
*
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2021-05-15    张众      创建文件
****************************************************************************/
#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <math.h>
#include "Constant.h"
//#include "main.h"

extern double debris_data[123][8];   //全局碎片信息

double estimate_T(const std::vector<int>& R, std::vector<double>& T, double initT);
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf);
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, int N);
double Dv_ij(int i, int j, double ts, double tf);
double Dv_All(double* T, int* R, int n);
double Dv_All(const std::vector<double>& T, int* R, int n);
double Dv_All(const std::vector<double>& T, const std::vector<int>& R, std::vector<double>& dv_sequence);
double domega(int i);

#endif

