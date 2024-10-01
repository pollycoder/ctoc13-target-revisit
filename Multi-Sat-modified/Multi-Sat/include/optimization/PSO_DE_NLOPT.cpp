#include "PSO_DE_NLOPT.h"
#include <iomanip>
#include <iostream>
#include "random_threads.h"

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
	std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut,
	double OmegaMin, double OmegaMax, double C1Min, double C1Max, double C2Min, double C2Max, double Vmax)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;

	int ct;
	double Omega, C1, C2, GBVAL;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[num_variable];
	double** popvel = NULL;
	popvel = new double* [Np];
	for (int i = 0;i < Np;i++)
		popvel[i] = new double[num_variable];

	std::vector<std::vector<double>> pop(Np);
	for (int i = 0;i < Np;i++)
		pop[i].resize(num_variable);
	double** PBPOS = new double* [Np];
	for (int i = 0;i < Np;i++)
		PBPOS[i] = new double[num_variable];

	for (int j = 0;j < Np;j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;
	for (int i = 0;i < Np;i++)
	{
		for (int j = 0;j < num_variable;j++)
		{
			pop[i][j] = realRand(0.0, 1.0);
			popvel[i][j] = 0.0;
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
		for (int i = 0;i < Np;i++)
		{
			double val = ObjFun(pop[i], grad, f_data);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0;j < num_variable;j++)
					PBPOS[i][j] = pop[i][j];
			}
			if (val < GBVAL)
			{
				GBVAL = val;
				for (int j = 0;j < num_variable;j++)
					GBPOS[j] = pop[i][j];
			}
		}
		Omega = OmegaMax - (OmegaMax - OmegaMin) / ItMax * ct;
		C1 = -(C1Max - C1Min) * ct / ItMax + C1Max;
		C2 = (C2Max - C2Min) * ct / ItMax + C2Min;
		for (int i = 0;i < Np;i++)
		{
			for (int j = 0;j < num_variable;j++)
			{
				popvel[i][j] = Omega * popvel[i][j] + C1 * (PBPOS[i][j] - pop[i][j]) * realRand(0.0, 1.0)
					+ C2 * (GBPOS[j] - pop[i][j]) * realRand(0.0, 1.0);
				if (popvel[i][j] > Vmax) popvel[i][j] = Vmax;//
				else if (popvel[i][j] < -Vmax) popvel[i][j] = -Vmax;//
				pop[i][j] += popvel[i][j];
				if ((pop[i][j] > 1.0) || (pop[i][j] < 0.0))
				{
					pop[i][j] = realRand(0.0, 1.0);
					popvel[i][j] = 0.0;
				}
			}
		}
		for (int j = 0;j < num_variable;j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			std::cout << "No. of iteration=" << ct << std::endl;
			for (int i = 0;i < num_variable;i++)
				std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
			std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
		}
		ct = ct + 1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
	for (int i = 0;i < Np;i++)
		delete[] popvel[i];
	delete[] popvel;
	for (int i = 0;i < Np;i++)
		delete[] PBPOS[i];
	delete[] PBPOS;
}

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
	std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut, double CR)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;
	if (CR < 0.0)
		CR = (18.0 + num_variable) / 40.0;
	if (CR > 0.9)
		CR = 0.9;

	int ct, r0, r1, r2;
	double F, GBVAL, lambda;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[num_variable];

	std::vector<std::vector<double>> pop(Np);
	for (int i = 0;i < Np;i++)
		pop[i].resize(num_variable);
	std::vector<std::vector<double>> u(Np);
	for (int i = 0;i < Np;i++)
		u[i].resize(num_variable);

	for (int j = 0;j < Np;j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;

	for (int i = 0;i < Np;i++)
	{
		for (int j = 0;j < num_variable;j++)
		{
			pop[i][j] = realRand(0.0, 1.0);
		}
		PBVAL[i] = ObjFun(pop[i], grad, f_data);
		if (PBVAL[i] < GBVAL)
		{
			GBVAL = PBVAL[i];
			for (int j = 0;j < num_variable;j++)
				GBPOS[j] = pop[i][j];
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
		lambda = exp(1 - ItMax / (ItMax + 1 - ct));
		F = 0.5 * pow(2, lambda);

		for (int i = 0; i < Np; i++)
		{
			r0 = intRand(0, Np - 1) % Np;
			r1 = intRand(0, Np - 1) % Np;
			r2 = intRand(0, Np - 1) % Np;
			for (int j = 0;j < num_variable;j++)
			{
				if (realRand(0.0, 1.0) <= CR)
					u[i][j] = pop[r0][j] + F * (pop[r1][j] - pop[r2][j]);
				else
					u[i][j] = pop[i][j];
				if ((u[i][j] > 1.0) || (u[i][j] < 0.0))
				{
					u[i][j] = realRand(0.0, 1.0);
				}
			}
			double val = ObjFun(u[i], grad, f_data);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0;j < num_variable;j++)
					pop[i][j] = u[i][j];
			}
		}

		for (int i = 0; i < Np; i++)
		{
			if (PBVAL[i] < GBVAL)
			{
				GBVAL = PBVAL[i];
				for (int j = 0;j < num_variable;j++)
					GBPOS[j] = pop[i][j];
			}
		}

		for (int j = 0;j < num_variable;j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			std::cout << "No. of iteration=" << ct << std::endl;
			for (int i = 0;i < num_variable;i++)
				std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
			std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
		}
		ct = ct + 1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
}

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
	std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut,
	double OmegaMin, double OmegaMax, double C1Min, double C1Max, double C2Min, double C2Max, double Vmax)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;

	int ct;
	double Omega, C1, C2, GBVAL;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[num_variable];
	double** popvel = NULL;
	popvel = new double* [Np];
	for (int i = 0;i < Np;i++)
		popvel[i] = new double[num_variable];

	std::vector<std::vector<double>> pop(Np);
	for (int i = 0;i < Np;i++)
		pop[i].resize(num_variable);
	double** PBPOS = new double* [Np];
	for (int i = 0;i < Np;i++)
		PBPOS[i] = new double[num_variable];

	for (int j = 0;j < Np;j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;
	for (int i = 0;i < Np;i++)
	{
		for (int j = 0;j < num_variable;j++)
		{
			pop[i][j] = realRand(0.0, 1.0);
			popvel[i][j] = 0.0;
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
#pragma omp parallel for schedule(dynamic)
		for (int i = 0;i < Np;i++)
		{
			double val = ObjFun(pop[i], grad, f_data);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0;j < num_variable;j++)
					PBPOS[i][j] = pop[i][j];
			}
#pragma omp critical 
			if (val < GBVAL)
			{
				GBVAL = val;
				for (int j = 0;j < num_variable;j++)
					GBPOS[j] = pop[i][j];
			}
		}
		Omega = OmegaMax - (OmegaMax - OmegaMin) / ItMax * ct;
		C1 = -(C1Max - C1Min) * ct / ItMax + C1Max;
		C2 = (C2Max - C2Min) * ct / ItMax + C2Min;
		for (int i = 0;i < Np;i++)
		{
			for (int j = 0;j < num_variable;j++)
			{
				popvel[i][j] = Omega * popvel[i][j] + C1 * (PBPOS[i][j] - pop[i][j]) * realRand(0.0, 1.0)
					+ C2 * (GBPOS[j] - pop[i][j]) * realRand(0.0, 1.0);
				if (popvel[i][j] > Vmax) popvel[i][j] = Vmax;//
				else if (popvel[i][j] < -Vmax) popvel[i][j] = -Vmax;//
				pop[i][j] += popvel[i][j];
				if ((pop[i][j] > 1.0) || (pop[i][j] < 0.0))
				{
					pop[i][j] = realRand(0.0, 1.0);
					popvel[i][j] = 0.0;
				}
			}
		}
		for (int j = 0;j < num_variable;j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			std::cout << "No. of iteration=" << ct << std::endl;
			for (int i = 0;i < num_variable;i++)
				std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
			std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
		}
		ct = ct + 1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
	for (int i = 0;i < Np;i++)
		delete[] popvel[i];
	delete[] popvel;
	for (int i = 0;i < Np;i++)
		delete[] PBPOS[i];
	delete[] PBPOS;
}

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
	std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut, double CR)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;
	if (CR < 0.0)
		CR = (18.0 + num_variable) / 40.0;
	if (CR > 0.9)
		CR = 0.9;

	int ct, r0, r1, r2;
	double F, GBVAL, lambda;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[num_variable];

	std::vector<std::vector<double>> pop(Np);
	for (int i = 0;i < Np;i++)
		pop[i].resize(num_variable);
	std::vector<std::vector<double>> u(Np);
	for (int i = 0;i < Np;i++)
		u[i].resize(num_variable);

	for (int j = 0;j < Np;j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0;i < Np;i++)
	{
		for (int j = 0;j < num_variable;j++)
		{
			pop[i][j] = realRand(0.0, 1.0);
		}
		PBVAL[i] = ObjFun(pop[i], grad, f_data);
#pragma omp critical
		if (PBVAL[i] < GBVAL)
		{
			GBVAL = PBVAL[i];
			for (int j = 0;j < num_variable;j++)
				GBPOS[j] = pop[i][j];
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
		lambda = exp(1 - ItMax / (ItMax + 1 - ct));
		F = 0.5 * pow(2, lambda);

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < Np; i++)
		{
			r0 = intRand(0, Np - 1) % Np;
			r1 = intRand(0, Np - 1) % Np;
			r2 = intRand(0, Np - 1) % Np;
			for (int j = 0;j < num_variable;j++)
			{
				if (realRand(0.0, 1.0) <= CR)
					u[i][j] = pop[r0][j] + F * (pop[r1][j] - pop[r2][j]);
				else
					u[i][j] = pop[i][j];
				if ((u[i][j] > 1.0) || (u[i][j] < 0.0))
				{
					u[i][j] = realRand(0.0, 1.0);
				}
			}
			double val = ObjFun(u[i], grad, f_data);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0;j < num_variable;j++)
					pop[i][j] = u[i][j];
			}
		}

		for (int i = 0; i < Np; i++)
		{
			if (PBVAL[i] < GBVAL)
			{
				GBVAL = PBVAL[i];
				for (int j = 0;j < num_variable;j++)
					GBPOS[j] = pop[i][j];
			}
		}

		for (int j = 0;j < num_variable;j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			std::cout << "No. of iteration=" << ct << std::endl;
			for (int i = 0;i < num_variable;i++)
				std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
			std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
		}
		ct = ct + 1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
}


void nlopt_main(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
	void* f_data, std::vector<double>& X, double& f, int num_variable, int Np, int maxstep )
{
	nlopt::opt opter(nlopt::LN_SBPLX, num_variable);      //����һ���Ż�����ʹ��LN_BOBYQA�㷨�����������Լ����ʹ�ö����Ż����ٶȸ���
											  //ʹ��LN_COBYLA�㷨(����Ҫ�ݶ����ܴ��������Լ��)��
											 // GN_DIRECT_L
											//����Ч������ LN_SBPLX�� GN_CRS2_LM  �ȽϺ�
	//int ineq_num = 3;  //������Լ������ά��
	//
	//std::vector<double> tol_ineq(ineq_num);
	//for (int i = 0; i < tol_ineq.size(); i++)  tol_ineq[i] = 1.0e-8;
	//opter.add_equality_mconstraint(inequality_constraint, f_data, tol_ineq);

	
	double tol = 1e-8;
	opter.set_xtol_abs(1e-8);
	opter.set_ftol_abs(1e-8);
	opter.set_force_stop(1e-8);
	opter.set_maxeval(maxstep);        //�Ż��ô�����ֹͣ

	std::vector<double> dx{ 0.003,0.0005 };
	std::vector<double> lb(X.size()), ub(X.size());     //�½� //�Ͻ�
	for (int i =0; i < X.size(); i++)
	{
		lb[i] = 0.0; ub[i] = 1.0;
	}
	opter.set_lower_bounds(lb);
	opter.set_upper_bounds(ub);
	//opter.set_initial_step(dx); //���ó�ʼ����,Ҳ���Բ�����
	opter.set_min_objective(ObjFun, f_data);   //ָ��



	
	nlopt::result res = opter.optimize(X, f);
	//int numers = opter.get_numevals();

}