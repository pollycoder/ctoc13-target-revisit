#include "PSO_DE_NLOPT.h"
#include <iomanip>
#include <iostream>
#include "random_threads.h"

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
	nlopt::opt opter(nlopt::LN_SBPLX, num_variable);      //定义一个优化器，使用LN_BOBYQA算法不处理非线性约束，使用二次优化，速度更快
											  //使用LN_COBYLA算法(不需要梯度且能处理非线性约束)，
											 // GN_DIRECT_L
											//测试效果来看 LN_SBPLX、 GN_CRS2_LM  比较好
	//int ineq_num = 3;  //非线性约束几个维度
	//
	//std::vector<double> tol_ineq(ineq_num);
	//for (int i = 0; i < tol_ineq.size(); i++)  tol_ineq[i] = 1.0e-8;
	//opter.add_equality_mconstraint(inequality_constraint, f_data, tol_ineq);

	
	double tol = 1e-8;
	opter.set_xtol_abs(1e-8);
	opter.set_ftol_abs(1e-8);
	opter.set_force_stop(1e-8);
	opter.set_maxeval(maxstep);        //优化该次数后停止

	std::vector<double> dx{ 0.003,0.0005 };
	std::vector<double> lb(X.size()), ub(X.size());     //下界 //上界
	for (int i =0; i < X.size(); i++)
	{
		lb[i] = 0.0; ub[i] = 1.0;
	}
	opter.set_lower_bounds(lb);
	opter.set_upper_bounds(ub);
	//opter.set_initial_step(dx); //设置初始步长,也可以不设置
	opter.set_min_objective(ObjFun, f_data);   //指标



	
	nlopt::result res = opter.optimize(X, f);
	//int numers = opter.get_numevals();

}