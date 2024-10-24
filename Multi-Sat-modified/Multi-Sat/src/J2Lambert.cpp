#include "J2Lambert.h"

#include "OrbitFun.h"

int shooting_func_lambert(int n, const double* x, double* fevc, int iflag, const double* para) {
    const double R1[3] = { para[0], para[1], para[2] };
    const double R2[3] = { para[3], para[4], para[5] };
    const double tof = para[6];
    const double v1[3] = { x[0], x[1], x[2] };
    
    double RV0[6] = { R1[0], R1[1], R1[2], v1[0], v1[1], v1[2]};
    double RVf[6];
    propagate_j2(RV0, RVf, 0.0, tof);
    
    fevc[0] = R2[0] - RVf[0];
    fevc[1] = R2[1] - RVf[1];
    fevc[2] = R2[2] - RVf[2];

    return 1;
}


int solve_lambert(const double* R1, const double* R2, const double& tof, double* v1vec, double* v2vec) {
    const double para[7] = { R1[0], R1[1], R1[2], R2[0], R2[1], R2[2], tof };
    double x[3] = { v1vec[0], v1vec[1], v1vec[2] };
    double fevc[3];
    MSolver A(3, 0);
    A.setMaxfevno(20);
    int flag = A.Solve(shooting_func_lambert, x, fevc, para);

    double err = V_Norm2(fevc, 3);
    if (err > 1e-3) return 0;

    for (int i = 0; i < 3; i++) v1vec[i] = x[i];
    double rv0[6] = { R1[0], R1[1], R1[2], v1vec[0], v1vec[1], v1vec[2] };
    double rvf[6];
    propagate_j2(rv0, rvf, 0.0, tof);

    for (int i = 0; i < 3; i++) v2vec[i] = rvf[i + 3];
    return 1;
}


void J2Lambert_short(int& flag, double* v1vec, double* v2vec, double& a, double& e, const double* R1, const double* R2, const double& TOF,
    const double& mu, int way=0, int N=0, int branch=0, int Maxiter=60, double tol=1e-11) {
    flag = lambert(v1vec, v2vec, a, e, R1, R2, TOF, mu, way, N, branch, Maxiter, tol);
    //if (flag != 1) std::cout << "Two-body Lambert failed." << std::endl;
    flag = solve_lambert(R1, R2, TOF, v1vec, v2vec);
    //if (flag == 0) std::cout << "J2 Lambert failed." << std::endl;
}

void J2Lambert_short(int&flag, double* RV1, double* RVf, const double* RV0, const double& TOF, const double& mu) {
    flag = 0;
    int flag_temp = 0;
    double R0[3] = { RV0[0], RV0[1], RV0[2] };
    double Rf[3] = { RVf[0], RVf[1], RVf[2] };
    double V0[3] = { RV0[3], RV0[4], RV0[5] };
    
    double a, e, v1temp[3], v2temp[3], dv1temp[3];
    double dvtemp = 1.0e10;
    double dv = 1.0e10;

    double radius = V_Norm2(R0, 3);
    double T = D2PI * sqrt(radius * radius * radius / mu_km_s);
    int Nmin = TOF / T - 1;
    int Nmax = TOF / T + 1;

    int way_opt, branch_opt;
    for (int N = Nmin; N <= Nmax; N++) {
        if (N < 0) continue;
        for (int way = 0; way < 2; way++) {
            for (int branch = 0; branch < 2; branch++) {
                J2Lambert_short(flag_temp, v1temp, v2temp, a, e, R0, Rf, TOF, mu, way, N, branch);
                if (flag_temp != 1) {
                    //std::cout << "J2求解失败" << std::endl;
                    continue;
                }
                if (a * (1 - e) - Re_km < 220.0) {
                    //std::cout << "近地点高度过低，peri = " << a * (1-e) - Re_km << " km" << std::endl;
                    continue;
                }
                if (a * (1 + e) - Re_km > 980.0) {
                    //std::cout << "远地点高度过高，apo = " << a * (1+e) - Re_km  << " km" << std::endl;
                    continue;
                }

                V_Minus(dv1temp, v1temp, V0, 3);
                dvtemp = V_Norm2(dv1temp, 3);
                if (dvtemp < dv) {
                    dv = dvtemp;
                    for (int i = 0; i < 3; i++) {
                        RVf[i + 3] = v2temp[i];
                        RV1[i] = RV0[i];
                        RV1[i + 3] = v1temp[i];
                    }
                    flag = flag_temp;
                    way_opt = way;
                    branch_opt = branch;
                }
            }
        }
    }

    //if (flag != 1) std::cout << "J2 Lambert failed." << std::endl;
    //else std::cout << "way = " << way_opt << " " << "branch = " << branch_opt << std::endl;*/
}

void test_lambert() {
    const double coe0[6] = { 7000.0, 0.0000, 0.0, 0.0, 0.0, 0.0 };
    const double coef[6] = { 7000.0, 0.5000, 0.3, 1.0, 2.0, 2.2 };
    double rv0[6], rv1[6], rvf[6], v1[3], v2[3], a, e; int flag;
    coe2rv(flag, rv0, coe0, mu_km_s);
    coe2rv(flag, rvf, coef, mu_km_s);

    double r1[3] = { rv0[0], rv0[1], rv0[2] };
    double r2[3] = { rvf[0], rvf[1], rvf[2] };

    J2Lambert_short(flag, v1, v2, a, e, r1, r2, 6000.0, mu_km_s, 0, 0, 0, 60, 1e-11);
    std::cout << "a = " << a << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "v1 = " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
    std::cout << "v2 = " << v2[0] << " " << v2[1] << " " << v2[2] << std::endl;

    J2Lambert_short(flag, rv1, rvf, rv0, 200000.0, mu_km_s);
    std::cout << "rv0 = " << rv0[0] << " " << rv0[1] << " " << rv0[2] << " " << rv0[3] << " " << rv0[4] << " " << rv0[5] << std::endl;
    std::cout << "rv1 = " << rv1[0] << " " << rv1[1] << " " << rv1[2] << " " << rv1[3] << " " << rv1[4] << " " << rv1[5] << std::endl;
    std::cout << "dv = " << rv1[3] - rv0[3] << " " << rv1[4] - rv0[4] << " " << rv1[5] - rv0[5] << std::endl;
}