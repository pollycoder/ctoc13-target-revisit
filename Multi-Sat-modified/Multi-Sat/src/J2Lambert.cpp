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


void J2Lambert_short(double* v1vec, double* v2vec, double& a, double& e, const double* R1, const double* R2, const double& TOF,
    const double& mu, int way=0, int N=0, int branch=0, int Maxiter=60, double tol=1e-11) {
    int flag = lambert(v1vec, v2vec, a, e, R1, R2, TOF, mu, way, N, branch);
    if (flag != 1) std::cout << "Two-body Lambert failed." << std::endl;
    flag = solve_lambert(R1, R2, TOF, v1vec, v2vec);
    if (flag == 0) std::cout << "J2 Lambert failed." << std::endl;
}

void test_lambert() {
    const double coe0[6] = { 7000.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const double coef[6] = { 7000.0, 0.0000, 0.0, 0.0, 0.0, 2.2 };
    double rv0[6], rvf[6], v1[3], v2[3], a, e; int flag;
    coe2rv(flag, rv0, coe0, mu_km_s);
    coe2rv(flag, rvf, coef, mu_km_s);

    double r1[3] = { rv0[0], rv0[1], rv0[2] };
    double r2[3] = { rvf[0], rvf[1], rvf[2] };

    J2Lambert_short(v1, v2, a, e, r1, r2, 6000.0, mu_km_s, 0, 0, 0, 60, 1e-11);
    std::cout << "a = " << a << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "v1 = " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
    std::cout << "v2 = " << v2[0] << " " << v2[1] << " " << v2[2] << std::endl;
}