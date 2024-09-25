#include"RungeKutta.h"
#include<stdlib.h>
namespace RK
{
double ArrayMax(const double* array, int dim)
{
	double temp=array[0];
	for(int i=0;i<dim;i++)
		temp=(temp<array[i])?array[i]:temp;
	return temp;
}

double ArrayMin(const double* array, int dim)
{
	double temp=array[0];
	for(int i=0;i<dim;i++)
		temp=(temp>array[i])?array[i]:temp;
	return temp;
}

double Max (double x, double y) 
{ return (x>y)?x:y; }

double Min (double x, double y) 
{ return (x<y)?x:y; }

double CopySign (double x, double y)
{ return (y < 0.0) ? ((x < 0.0) ? x : -x) : ((x > 0.0) ? x : -x); }

int CopySign (int x, int y)
{ return (y < 0.0) ? ((x < 0.0) ? x : -x) : ((x > 0.0) ? x : -x); }

//RKF78误差估计
double RKF78ErrorTerm(int k, int dim, double* RKFWork)
{
	//RKF78的系数
	static const double RKFe[13]=
		{ 41.0/840.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 41.0/840.0, -41.0/840.0, -41.0/840.0 };
    double sum = 0.0;
    for (int i = 0; i <= 12; i++)
        sum += RKFe[i] * RKFWork[i*dim+k];
    return sum;
}

//RKF78右函数值计算
void RKF78Step(void (*Model)(double t, const double* y, double* yp, const double* para), const double* y, const double* para,
			 double t, double h, double* s, int dim, double* RKFWork)
{
	//RKF78的系数
	static const double RKFa[13]=
		{ 0.0, 2.0/27.0,  1.0/9.0,  1.0/6.0,  5.0/12.0, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0,  0.0, 1.0};
	static const double RKFb[13][12]=
		{
			{  0.0,            0.0,       0.0,         0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  2.0/27.0,       0.0,       0.0,         0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  1.0/36.0,       1.0/12.0,  0.0,         0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  1.0/24.0,       0.0,       1.0/8.0,     0.0,           0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  5.0/12.0,       0.0,       -25.0/16.0,  25.0/16.0,     0.0,           0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  1.0/20.0,       0.0,       0.0,         1.0/4.0,       1.0/5.0,       0.0,         0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  -25.0/108.0,    0.0,       0.0,         125.0/108.0,   -65.0/27.0,    125.0/54.0,  0.0,           0.0,       0.0,        0.0,       0.0, 0.0, },
			{  31.0/300.0,     0.0,       0.0,         0.0,           61.0/225.0,    -2.0/9.0,    13.0/900.0,    0.0,       0.0,        0.0,       0.0, 0.0, },
			{  2.0,            0.0,       0.0,         -53.0/6.0,     704.0/45.0,    -107.0/9.0,  67.0/90.0,     3.0,       0.0,        0.0,       0.0, 0.0, },
			{  -91.0/108.0,    0.0,       0.0,         23.0/108.0,    -976.0/135.0,  311.0/54.0,  -19.0/60.0,    17.0/6.0,  -1.0/12.0,  0.0,       0.0, 0.0, },
			{  2383.0/4100.0,  0.0,       0.0,         -341.0/164.0,  4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0, 0.0, 0.0, },
			{  3.0/205.0,      0.0,       0.0,         0.0,           0.0,           -6.0/41.0,   -3.0/205.0,    -3.0/41.0, 3.0/41.0,   6.0/41.0,  0.0, 0.0, },
			{  -1777.0/4100.0, 0.0,       0.0,         -341.0/164.0,  4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0.0, 1.0, }
	   };
	static const double RKFc[13]=
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0};
	
    int hi=dim-1;	
	for(int j=1;j<=12;j++)
	{
        for(int i =0;i<=hi;i++)
		{
            double x = 0.0;
            for (int m = 0; m < j; m++) 
                x += RKFb[j][m] * RKFWork[m*dim+i];
            s[i] = x * h + y[i];
        }
        Model(t + RKFa[j] * h, s, &RKFWork[j*dim],para);
    }

    for (int i = 0; i <= hi; i++)
	{
        double x = 0.0;
        for (int j = 0; j <= 12; j++) 
            x += RKFc[j] * RKFWork[j*dim+i];
        s[i] = h * x + y[i];
    }	
}

//RKF78积分函数,t以初值y积分至tout,结果存于y中.
int RKF78(void (*Model)(double t, const double* y, double* yp, const double* para), double* y, const double* para,
		  double& t, double tout, int dim, double RelTol, double* AbsTol, int& RKFiflag, double* RKFWork,
		  int RKFmaxnfe, FILE* fid)
{	
    // Get machine epsilon
 
    const double eps = 2.2204460492503131e-016,
                 u26 = 26*eps;               
    
    const double remin = 1e-15;	

    int mflag = abs(RKFiflag);
	int RKFnfe, RKFkop, RKFinit=0, RKFkflag=0, RKFjflag=0;
	double RKFh, RKFsavre=0.0, RKFsavae=0.0;
	
	if(fabs(t - tout)<1.0E-14)
	{
		RKFiflag*=2;
		t=tout;
		if(fid)
		{
			fprintf(fid,"%22.14e",t);
			for(int i=0;i<dim;i++)
				fprintf(fid,"%24.14e",y[i]);
			fprintf(fid,"\n");
			fflush(fid);
		}
		return 1;
	}

    if (dim < 1 || RelTol < 0.0 || ArrayMin(AbsTol, dim) < 0.0 ||  mflag == 0  
        || mflag > 8 ) 
	{
        RKFiflag = 8;
        return 0;
    }

    double dt,rer,scale,hmin,eeoet,ee,et,esttol,s,ae,tol = 0;
	double iorder=8.0;
    int output, hfaild, k;
	int NumNode=0;
/*	if(fid)
	{
		NumNode++;
		fprintf(fid,"%10d%10d\n",0,dim+1);
		fprintf(fid,"%22.14e",t);
		for(int i=0;i<dim;i++)
			fprintf(fid,"%24.14e",y[i]);
		fprintf(fid,"\n");
		fflush(fid);
	}*/
    int lo = 0, hi = dim-1;    
    int gflag = 0;
         
    if (RKFiflag == 3||(mflag == 2 && (RKFinit == 0 || RKFkflag == 2))) 
	{
        gflag = 1;
        goto next;
    }

    if (RKFiflag == 4 || (RKFkflag == 4 && mflag == 2))
	{
        RKFnfe = 0;
        if (mflag != 2) gflag = 1;
        goto next;
    }
        
    if ((RKFkflag == 5 && ArrayMin(AbsTol,dim) == 0.0)
        || (RKFkflag == 6 && RelTol < RKFsavre && ArrayMin(AbsTol,dim) < RKFsavae))
	{
        RKFiflag = 9;		
        goto final;
    }
   
  next:

    if (gflag) 
	{
        RKFiflag = RKFjflag;
        if (RKFkflag == 3) mflag = abs(RKFiflag);
    }

    RKFjflag = RKFiflag;
    RKFkflag = 0;
    RKFsavre = RelTol;
    RKFsavae = ArrayMin(AbsTol,dim);   
    
    rer = 2 * eps + remin;    
    if (RelTol < rer) 
	{
    //    RelTol = rer;
        RKFiflag = RKFkflag = 3;
        goto final;
    }
    
    gflag = 0;
    dt = tout - t;

    if (mflag == 1)
	{
        RKFinit = 0;
        RKFkop = 0;
        gflag = 1;
        Model(t,y,&RKFWork[0],para);  // call function
        RKFnfe = 1;
        if (fabs(t - tout)<1.0E-14) 
		{
            RKFiflag = 2;
            goto final;
        }
    }

    if (RKFinit == 0 || gflag)
	{
        RKFinit = 1;
        RKFh = fabs(dt);
        double ypk;
        for (int k = lo; k <= hi; k++)
		{
            tol = Max(RelTol*fabs(y[k]), AbsTol[k]);//RelTol * fabs(y(k)) + AbsTol;//
            if (tol > 0.0) 
			{
                ypk = fabs(RKFWork[k]);
				double RKFh8pow=pow(RKFh,iorder);//RKFh*RKFh;
		//		RKFh8pow*=(RKFh8pow*=RKFh8pow);
                if (ypk * RKFh8pow > tol)
				{
					double temp=tol/ypk;					
                    RKFh = pow(temp,1.0/iorder);//sqrt(sqrt(sqrt(temp)));
				}
            }
        }

        if (tol <= 0.0) RKFh = 0.0;
        ypk = Max(fabs(dt),fabs(t));
        RKFh = Max(RKFh, u26 * ypk);
        RKFjflag = CopySign(2,RKFiflag);
    }

    // Set stepsize for integration in the direction from t to tout

    RKFh = CopySign(RKFh,dt);

    // Test to see if this routine is being severely impacted by too many
    // output points

    if (fabs(RKFh) >= 2*fabs(dt)) RKFkop++;

    if (RKFkop == 100) {
        RKFkop = 0;
        RKFiflag = 7;
        goto final;
    }

    if (fabs(dt) <= u26 * fabs(t)) 
	{
        // If too close to output point,extrapolate and return
        for (int k = lo; k <= hi; k++)
            y[k] += dt * RKFWork[k];

        Model(tout,y,&RKFWork[0],para);
        RKFnfe++;
        t = tout;
        RKFiflag = 2;
        goto final;
    }

    // Initialize output point indicator

    output = false;

    // To avoid premature underflow in the error tolerance function,
    // scale the error tolerances

    scale = 2.0 / RelTol;
    ae = scale * ArrayMin(AbsTol,dim); 

    // Step by step integration - as an endless loop over steps

    for (;;) 
	{ 

        hfaild = 0;

        // Set smallest allowable stepsize

        hmin = u26 * fabs(t);
        dt = tout - t;
        if (fabs(dt) < 2.0 * fabs(RKFh))
		{
            if (fabs(dt) <= fabs(RKFh)) 
			{

                // The next successful step will complete the 
                // integration to the output point

                output = true;
                RKFh = dt;
            }
			else
                RKFh = 0.5 * dt;
        }
        
        if (RKFnfe > RKFmaxnfe)
		{
            RKFiflag = RKFkflag = 4;
            goto final;
        }
        
    step:

        RKF78Step(Model,y, para,t, RKFh, &RKFWork[13*dim], dim, RKFWork);      
        for (int i = lo; i <= hi; i++) RKFWork[dim+i] = RKFWork[13*dim+i];
        RKFnfe += 8;
        eeoet = 0.0;
        for (k = lo; k <= hi; k++)
		{

            et = fabs(y[k]) + fabs(RKFWork[dim+k]) + ae;
            // Inappropriate error tolerance

            if (et <= 0.0) 
			{
                RKFiflag = 5;
                goto final;
            }
            
            ee = fabs( RKF78ErrorTerm(k,dim,RKFWork));
//			ee=Max(ee,pow(RKFh,8.0));
            eeoet = Max(eeoet,ee/et);
        }
	//	eeoet=Max(eeoet,1e-14);
//		printf("%.6e,%.6e\n",eeoet,RKFh);
        
        esttol = fabs(RKFh) * eeoet * scale;
        
        if (esttol > 1.0)
		{
           
            hfaild = true;
            output = false;
            s = 0.1;
            if (esttol < pow(9.0,iorder)) s = 0.9 / pow(esttol, 1.0/iorder);//sqrt(sqrt(sqrt(esttol)));//43046721.0
            RKFh *= s;
            if (fabs(RKFh) > hmin) goto step; // loop

            // Requested error unattainable at smallest allowable stepsize

            RKFiflag = RKFkflag = 6;
            goto final;
        }

        // Successful step
        // Store solution at t+h and evaluate derivatives there

        t += RKFh;
//		cout<<t<<","<<RKFh<<endl;
        for (k = lo; k <= hi; k++) y[k] = RKFWork[dim+k];

        Model(t,y,&RKFWork[0],para);
        RKFnfe++;     
        s = 5.0;
        if (esttol > pow(0.9/iorder,iorder)) s = 0.9 / pow(esttol, 1.0/iorder);//sqrt(sqrt(sqrt(esttol)));//2.565784513950347900390625E-8
        if (hfaild) s = Min(1.0, s);
        RKFh = CopySign(Max(hmin,s*fabs(RKFh)), RKFh);   
		
		if(fid)
		{
			NumNode++;
			fprintf(fid,"%22.14e",t);
			for(int i=0;i<dim;i++)
				fprintf(fid,"%24.14e",y[i]);
			fprintf(fid,"\n");
			fflush(fid);
		}
        

        if (output) 
		{
            t = tout;
			RKFiflag = 2;		
			goto final;
        }

        if (RKFiflag <= 0)
		{ // one-step mode
            RKFiflag = -2;
            goto final;
        }
    }
final:
/*	if(fid)
	{
		rewind(fid);
		fprintf(fid,"%10d",NumNode);
	}*/
			
//	for(i=0;i<14;i++)
//		delete[] RKFWork[i];
//	delete[] RKFWork;
	return NumNode;
}

//wa:2*dim
void RK4(void (*Model)(double t, const double* y, double* yp, const double* para),
			 const double* y, const double* para, double t, double h, double* s, int dim, double* wa)
{
	int i;
	if(fabs(h)<=0.0)
	{
		for(i=0;i<dim;i++)
			s[i]=y[i];
		return;
	}
//	double* k1=new double[dim];
//	double* k2=new double[dim];
	Model(t,y,&wa[dim],para);
	for(i=0;i<dim;i++)
	{
		s[i]=h*wa[dim+i];
		wa[i]=0.5*h*wa[dim+i]+y[i];
	}
	Model(t+0.5*h,wa,&wa[dim],para);
	for(i=0;i<dim;i++)
	{
		s[i]+=2.0*h*wa[dim+i];
		wa[i]=0.5*h*wa[dim+i]+y[i];
	}
	Model(t+0.5*h,wa,&wa[dim],para);
	for(i=0;i<dim;i++)
	{
		s[i]+=2.0*h*wa[dim+i];
		wa[i]=h*wa[dim+i]+y[i];
	}
	Model(t+h,wa,&wa[dim],para);
	for(i=0;i<dim;i++)
		s[i]=(s[i]+h*wa[dim+i])/6.0+y[i];

//	delete[] k1;
//	delete[] k2;
}

int RK4(void (*Model)(double t, const double* y, double* yp, const double* para),
			 double* y, const double* para, double& tin, double tout, double h, int dim, double* wa, FILE* fid)
{
	double dt;
	int Num=0,i;
	bool final=false;
	if(tout<tin) h=-h;
	while(!final)
	{
		dt=h;
		if(fabs(h)>=fabs(tin-tout))
		{
			dt=tout-tin;
			final=true;
		}
		RK4(Model, y, para, tin, dt, &wa[2*dim], dim, wa);
		for(i=0;i<dim;i++) y[i]=wa[2*dim+i];
		tin+=dt;
		if(fid)
		{
			Num++;
			fprintf(fid,"%22.14e",tin);
			for(int i=0;i<dim;i++)
				fprintf(fid,"%24.14e",y[i]);
			fprintf(fid,"\n");
			fflush(fid);
		}		
	}
	return Num;
}


//RKF45误差估计
double RKF45ErrorTerm(int k, int dim, double* RKFWork)
{
	//RKF78的系数
	static const double RKFe[6]=
		{-1.0/360.0,  0.0, 128.0/4275.0, 2197.0/75240.0, -1.0/50.0, -2.0/55.0};
    double sum = 0.0;
    for (int i = 0; i <= 5; i++)
        sum += RKFe[i] * RKFWork[i*dim+k];
    return sum;
}

//RKF45右函数值计算
void RKF45Step(void (*Model)(double t, const double* y, double* yp, const double* para), const double* y, const double* para,
			 double t, double h, double* s, int dim, double* RKFWork)
{
	//RKF78的系数
	static const double RKFa[6]=
		{ 0.0, 1.0/4.0,  3.0/8.0,  12.0/13.0,  1.0,  1.0/2.0};
	static const double RKFb[6][5]=
		{
	{  0.0,            0.0,           0.0,            0.0,            0.0       },
    {  1.0/4.0,        0.0,           0.0,            0.0,            0.0       },
    {  3.0/32.0,       9.0/32.0,      0.0,            0.0,            0.0       },
    {  1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,  0.0,            0.0       },
    {  439.0/216.0,   -8.0,           3680.0/513.0,  -845.0/4104.0,   0.0       },
    { -8.0/27.0,       2.0,          -3544.0/2565.0,  1859.0/4104.0, -11.0/40.0 }
	   };
	static const double RKFc[6]=
		{ 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
	
    int hi=dim-1;	
	for(int j=1;j<=5;j++)
	{
        for(int i =0;i<=hi;i++)
		{
            double x = 0.0;
            for (int m = 0; m < j; m++) 
                x += RKFb[j][m] * RKFWork[m*dim+i];
            s[i] = x * h + y[i];
        }
        Model(t + RKFa[j] * h, s, &RKFWork[j*dim],para);
    }

    for (int i = 0; i <= hi; i++)
	{
        double x = 0.0;
        for (int j = 0; j <= 5; j++) 
            x += RKFc[j] * RKFWork[j*dim+i];
        s[i] = h * x + y[i];
    }	
}

//RKF78积分函数,t以初值y积分至tout,结果存于y中.
int RKF45(void (*Model)(double t, const double* y, double* yp, const double* para), double* y, const double* para,
		  double& t, double tout, int dim, double RelTol, double* AbsTol, int& RKFiflag, double* RKFWork,
		  int RKFmaxnfe, FILE* fid)
{	
    // Get machine epsilon
 
    const double eps = 2.2204460492503131e-016,
                 u26 = 26*eps;               
    
    const double remin = 1e-15;	

    int mflag = abs(RKFiflag);
	int RKFnfe, RKFkop, RKFinit=0, RKFkflag=0, RKFjflag=0;
	double RKFh, RKFsavre=0.0, RKFsavae=0.0;

	if(fabs(t - tout)<1.0E-14)
	{
		RKFiflag*=2;
		t=tout;
		if(fid)
		{
			fprintf(fid,"%22.14e",t);
			for(int i=0;i<dim;i++)
				fprintf(fid,"%24.14e",y[i]);
			fprintf(fid,"\n");
			fflush(fid);
		}
		return 1;
	}

    if (dim < 1 || RelTol < 0.0 || ArrayMin(AbsTol, dim) < 0.0 ||  mflag == 0  
        || mflag > 8 || ((fabs(t - tout)<1.0E-14) && (RKFkflag != 3))) 
	{
        RKFiflag = 8;
        return 0;
    }

    double dt,rer,scale,hmin,eeoet,ee,et,esttol,s,ae,tol = 0;
    int output, hfaild, k;
	int NumNode=0;

    int lo = 0, hi = dim-1;    
    int gflag = 0;
         
    if (RKFiflag == 3||(mflag == 2 && (RKFinit == 0 || RKFkflag == 2))) 
	{
        gflag = 1;
        goto next;
    }

    if (RKFiflag == 4 || (RKFkflag == 4 && mflag == 2))
	{
        RKFnfe = 0;
        if (mflag != 2) gflag = 1;
        goto next;
    }
        
    if ((RKFkflag == 5 && ArrayMin(AbsTol,dim) == 0.0)
        || (RKFkflag == 6 && RelTol < RKFsavre && ArrayMin(AbsTol,dim) < RKFsavae))
	{
        RKFiflag = 9;		
        goto final;
    }
   
  next:

    if (gflag) 
	{
        RKFiflag = RKFjflag;
        if (RKFkflag == 3) mflag = abs(RKFiflag);
    }

    RKFjflag = RKFiflag;
    RKFkflag = 0;
    RKFsavre = RelTol;
    RKFsavae = ArrayMin(AbsTol,dim);   
    
    rer = 2 * eps + remin;    
    if (RelTol < rer) 
	{
    //    RelTol = rer;
        RKFiflag = RKFkflag = 3;
        goto final;
    }
    
    gflag = 0;
    dt = tout - t;

    if (mflag == 1)
	{
        RKFinit = 0;
        RKFkop = 0;
        gflag = 1;
        Model(t,y,&RKFWork[0],para);  // call function
        RKFnfe = 1;
        if (fabs(t - tout)<1.0E-14) 
		{
            RKFiflag = 2;
            goto final;
        }
    }

    if (RKFinit == 0 || gflag)
	{
        RKFinit = 1;
        RKFh = fabs(dt);
        double ypk;
        for (int k = lo; k <= hi; k++)
		{
            tol = Max(RelTol*fabs(y[k]), AbsTol[k]);//RelTol * fabs(y(k)) + AbsTol;//
            if (tol > 0.0) 
			{
                ypk = fabs(RKFWork[k]);
				double RKFh5pow=RKFh*RKFh*RKFh*RKFh*RKFh;				
                if (ypk * RKFh5pow > tol)
				{
					double temp=tol/ypk;					
                    RKFh = pow(temp,0.2);
				}
            }
        }

        if (tol <= 0.0) RKFh = 0.0;
        ypk = Max(fabs(dt),fabs(t));
        RKFh = Max(RKFh, u26 * ypk);
        RKFjflag = CopySign(2,RKFiflag);
    }

    // Set stepsize for integration in the direction from t to tout

    RKFh = CopySign(RKFh,dt);

    // Test to see if this routine is being severely impacted by too many
    // output points

    if (fabs(RKFh) >= 2*fabs(dt)) RKFkop++;

    if (RKFkop == 100) {
        RKFkop = 0;
        RKFiflag = 7;
        goto final;
    }

    if (fabs(dt) <= u26 * fabs(t)) 
	{
        // If too close to output point,extrapolate and return
        for (int k = lo; k <= hi; k++)
            y[k] += dt * RKFWork[k];

        Model(tout,y,&RKFWork[0],para);
        RKFnfe++;
        t = tout;
        RKFiflag = 2;
        goto final;
    }

    // Initialize output point indicator

    output = false;

    // To avoid premature underflow in the error tolerance function,
    // scale the error tolerances

    scale = 2.0 / RelTol;
    ae = scale * ArrayMin(AbsTol,dim); 

    // Step by step integration - as an endless loop over steps

    for (;;) 
	{ 

        hfaild = 0;

        // Set smallest allowable stepsize

        hmin = u26 * fabs(t);
        dt = tout - t;
        if (fabs(dt) < 2.0 * fabs(RKFh))
		{
            if (fabs(dt) <= fabs(RKFh)) 
			{

                // The next successful step will complete the 
                // integration to the output point

                output = true;
                RKFh = dt;
            }
			else
                RKFh = 0.5 * dt;
        }
        
        if (RKFnfe > RKFmaxnfe)
		{
            RKFiflag = RKFkflag = 4;
            goto final;
        }
        
    step:

        RKF45Step(Model,y, para,t, RKFh, &RKFWork[6*dim], dim, RKFWork);      
        for (int i = lo; i <= hi; i++) RKFWork[dim+i] = RKFWork[6*dim+i];
        RKFnfe += 8;
        eeoet = 0.0;
        for (k = lo; k <= hi; k++)
		{

            et = fabs(y[k]) + fabs(RKFWork[dim+k]) + ae;
            // Inappropriate error tolerance

            if (et <= 0.0) 
			{
                RKFiflag = 5;
                goto final;
            }
            
            ee = fabs( RKF45ErrorTerm(k,dim,RKFWork));
			eeoet = Max(eeoet,ee/et);			
        }
//		eeoet=Max(eeoet,pow(RKFh,5.0));//后改，按上一式估计误差可能为零，会误导对步长的选取
//		printf("%.6e,%.6e\n",eeoet,RKFh);
        
        esttol = fabs(RKFh) * eeoet * scale;
        
        if (esttol > 1.0)
		{
           
            hfaild = true;
            output = false;
            s = 0.1;
            if (esttol < 5904.0) s = 0.9 / pow(esttol,0.2);//pow(esttol, 0.125);
            RKFh *= s;
            if (fabs(RKFh) > hmin) goto step; // loop

            // Requested error unattainable at smallest allowable stepsize

            RKFiflag = RKFkflag = 6;
            goto final;
        }

        // Successful step
        // Store solution at t+h and evaluate derivatives there

        t += RKFh;
//		cout<<t<<","<<RKFh<<endl;
        for (k = lo; k <= hi; k++) y[k] = RKFWork[dim+k];

        Model(t,y,&RKFWork[0],para);
        RKFnfe++;     
        s = 5.0;
        if (esttol > 1.889568E-4) s = 0.9 / pow(esttol,0.2);//pow(esttol, 0.125);
        if (hfaild) s = Min(1.0, s);
        RKFh = CopySign(Max(hmin,s*fabs(RKFh)), RKFh);   
		
		if(fid)
		{
			NumNode++;
			fprintf(fid,"%22.14e",t);
			for(int i=0;i<dim;i++)
				fprintf(fid,"%24.14e",y[i]);
			fprintf(fid,"\n");
			fflush(fid);
		}
        

        if (output) 
		{
            t = tout;
			RKFiflag = 2;		
			goto final;
        }

        if (RKFiflag <= 0)
		{ // one-step mode
            RKFiflag = -2;
            goto final;
        }
    }
final:
/*	if(fid)
	{
		rewind(fid);
		fprintf(fid,"%10d",NumNode);
	}*/
			
//	for(i=0;i<14;i++)
//		delete[] RKFWork[i];
//	delete[] RKFWork;
	return NumNode;
}
}
/*
void test(double t, const double* y, double* yp, const double* para)
{
	double s=cos(100.0*t);//cos(30.0*t);//1.0/(t+1.0);//
	yp[0]=s;//
//	if(s<0.0)
//		yp[0]=1.0;
//	yp[0]=y[2];
//	yp[1]=y[3];
//	double r=sqrt(y[0]*y[0]+y[1]*y[1]);
//	yp[2]=-3.986e14/(r*r*r)*y[0];
//	yp[3]=-3.986e14/(r*r*r)*y[1];
}
#include<time.h>
void main()
{
	clock_t start, stop;
	start=clock();
	double RelTol=1e-12,AbsTol[1]={1e-12};//,1e-3,1e-6,1e-6};	
	double wa[56],x[1]={0.0};//4e7,0.0,0.0,3.156738823532920e+003};
	double para[1],ta=0.0,tb=20.0;//7.961615652634385e+004;
	int Num,RKFiflag=1;
	FILE* fid=fopen("test1.txt","w");
	Num=RKSolver(test, x, para, ta, tb, 1, RelTol, AbsTol, RKFiflag, wa, 1000000,fid);
//	Num=RK::RK4(test,x, para, ta, tb, 1e-4, 1, wa, fid);
	printf("%.12f\n",x[0]);
	printf("%d\n",RKFiflag);
	fclose(fid);
	stop=clock();
	printf("%.12f\n",(double)(stop-start)/CLOCKS_PER_SEC);	
}*/