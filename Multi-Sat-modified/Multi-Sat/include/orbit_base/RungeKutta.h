#ifndef _RUNGEKUTTAF_H_
#define _RUNGEKUTTAF_H_
#include<iostream>
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"Constant.h"
using namespace std;

//Rugge-Kutta����
namespace RK{
//	ע�⣺ÿ�ε���RKF78���б䲽������ʱǰ���轫RKFiflag��Ϊ1�����ֳɹ�����ֵΪ2��RKF78���ص�����Ϊ�䲽�����ֵĲ���
//����˵����
	//Model:������������yp=f(t,y,para)������yΪdimά״̬��,ypΪ״̬����һ�׵���,paraΪ����yp�����һЩ�����������鱣��
	//y:dimά״̬���ĳ�ֵ�����ֳɹ��󣬱�������״ֵ̬
	//para:����yp�����һЩ�����������鱣��
	//tin:������ʼʱ�䣬���ֳɹ�����ֵ��Ϊ������ֹʱ��tout
	//tout:������ֹʱ��
	//dim:ά��
	//RelTol:������
	//AbsTol:��������״̬��ά��һ��
	//RKFiflag���룺1ȫ���俪ʼ���֣�-1��������ʼ����
	//		�����2ȫ������ֳɹ���-2���������ֳɹ���3������̫С��4�����������ޣ�5�������̫С��6�Ѵﵽ��С��������7�����Һ����Ĵ������ޣ�8������Ĳ�����9�������
	//RKFWork:һά��������14*dim��Ԫ�أ����û��������ȷ��䡣����״̬��Ϊ6ά������Ҫ���ִ�СΪ84��RKFWork��
	//RKFmaxnfe:���������������һ���ϴ���������ɣ���1000000��
	//fid:������ֹ��̵��ļ��������ΪNULL������ֵ���򲻱������ֵ����ֱ�Ӵ�y�еõ�����ֵ
int RKF78(void (*Model)(double t, const double* y, double* yp, const double* para), double* y, const double* para,
		  double& tin, double tout, int dim, double RelTol, double* AbsTol, int& RKFiflag, double* RKFWork,
		  int RKFmaxnfe=100000, FILE* fid=NULL);
//wa:2*dim
void RK4(void (*Model)(double t, const double* y, double* yp, const double* para),
			 const double* y, const double* para, double t, double h, double* s, int dim, double* wa);
//wa:3*dim
int RK4(void (*Model)(double t, const double* y, double* yp, const double* para),
			 double* y, const double* para, double& tin, double tout, double h, int dim, double* wa, FILE* fid=NULL);

//RKFWork:һά��������7Xdim��Ԫ�ء�RKFiflag,����RKF45ǰ����Ϊ1���Ա����ȫ����䲽�����֣����ֳɹ�����ֵΪ2
int RKF45(void (*Model)(double t, const double* y, double* yp, const double* para), double* y, const double* para,
		  double& tin, double tout, int dim, double RelTol, double* AbsTol, int& RKFiflag, double* RKFWork,
		  int RKFmaxnfe=100000, FILE* fid=NULL);
}

#endif

