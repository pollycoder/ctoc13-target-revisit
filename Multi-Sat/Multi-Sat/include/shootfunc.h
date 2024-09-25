/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: shootfunc.h
* 内容简述：打靶得到真实速度增量
*
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2020-10-04    张众      创建文件
****************************************************************************/

#ifndef _SHOOTFUNC_H
#define _SHOOTFUNC_H

#include "MinpackSolver.h"
#include "propagate.h"
#include "Constant.h"
#include "OrbitMath.h"
//#include "MyPoint.h"

using namespace MinpackSolver;

/****************************************************************************
* 函数名   : shoot_func()
* 功  能   : 打靶函数
* 输 入    : x: 0:dvx 1:dvy 2:dvz
* 参数     : para: 0-5:rv0 6:mass 7:if_drag(考虑为1.0 不考虑为0.0) 8:t0 9-11:地面点在tf下的地固系单位矢量 12:tf
* 输    出 : fevc: 0:0 1:0.0 2:0 3:0
****************************************************************************/
int shoot_func(int n, const double* x, double* fevc, int iflag, const double* para);

/****************************************************************************
* 函数名   : shoot_func_2()
* 功  能   : 打靶函数
* 输 入    : x: 0:dv
* 参数     : para: 0-5:rv0 6:mass 7:if_drag(考虑为1.0 不考虑为0.0) 8:t0 9-10:地面点经纬度（rad） 11:tf 12-14:dv的方向
* 输    出 : fevc: 0:经度差
****************************************************************************/
int shoot_func_2(int n, const double* x, double* fevc, int iflag, const double* para);

void acc_validate(const double* x, const double* para, double& t_access, double* rvf);
int acc_validate_total(const double& longitude, const double& latitude, const double* RV, double m0, double* dv, const double t0, double& tf, double* rvf, double& mf);

/****************************************************************************
* 函数名   : solve_shoot_dv()
* 功  能   : 通过打靶方程求解在大气下真实dv与tf
* 输入输出 :  RV为变轨前位置速度，m为没有变轨前的质量（输出为最后变轨后的质量），
              dv为初始预测的dv（输出为最后的dv），tf为初始预测的tf，不再改变
****************************************************************************/
//int solve_shoot_dv(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double& tf);

/****************************************************************************
* 函数名   : solve_trans_dv()
* 功  能   : 通过打靶方程求解在大气下真实dv与tf
* 输入输出 :  RV为变轨前位置速度，m为没有变轨前的质量（输出为最后变轨后的质量），
			  dv为初始预测的dv（输出为最后的dv），tf为初始预测的tf，改变
* 参数     : para: 0-5:rv0 6:mass 7:if_drag(考虑为1.0 不考虑为0.0) 8:t0 9-10:地面点经纬度（rad） 11:tf 12-14:dv的方向
****************************************************************************/
//int solve_trans_dv(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double tf);
int solve_trans_dv(const double& longitude, const double& latitude, const double* RV, double& m, double* dv, const double t0, double tf);

/****************************************************************************
* 函数名   : shoot_func_ToRSAT
* 功  能   : 打靶函数,到中继星
* 输 入    : x: 0:dvx 1:dvy 2:dvz
* 参数     : para: 0-5:rv0 6:mass 7:if_drag(考虑为1.0 不考虑为0.0) 8:t0 9-11:rvf 12:tf
* 输    出 : fevc: 0:0 1:0.0 2:0 3:0
****************************************************************************/
int shoot_func_ToRSAT(int n, const double* x, double* fevc, int iflag, const double* para);

/****************************************************************************
* 函数名   : solve_shoot_dv_ToRSAT()
* 功  能   : 通过打靶方程求解转移至中继星位置的准确的速度增量
* 输入输出 :  RV为变轨前位置速度，Rf为目标位置，m为没有变轨前的质量（输出为最后变轨后的质量），
			  dv为初始预测的dv（输出为最后的dv），tf为初始预测的tf，不再改变
****************************************************************************/
int solve_shoot_dv_ToRSAT(bool& flag, const double* RV, double* Rf, double& m, double* dv, const double t0, double& tf);

int solve_shoot_dv_ToRSAT2(bool& flag, const double* RV, double* Rf, double& m, double* dv, const double t0, double& tf);

int shoot_func3(int n, const double* x, double* fevc, int iflag, const double* para);

//int solve_shoot_dv3(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double& tf);

#endif