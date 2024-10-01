#include "visibility_4_targets.h"

#include <algorithm>

// 矩形视角的可见性
bool is_target_visible(const double* rv_sat, const double* rv_target, double half_cone_angle_rad) {
    const double r_sat[3] = { rv_sat[0], rv_sat[1], rv_sat[2] };
    const double v_sat[3] = { rv_sat[3], rv_sat[4], rv_sat[5] };
    const double r_dot_v_sat = V_Dot(r_sat, v_sat, 3);
    const double r_norm_sat = V_Norm2(r_sat, 3);
    const double v_norm_sat = V_Norm2(v_sat, 3);
    double eVec_r[3];
    V_Divid(eVec_r, r_sat, r_norm_sat, 3);//
    double v_sat_vert[3];
    V_Multi(v_sat_vert, eVec_r, r_dot_v_sat / r_norm_sat, 3);
    double v_sat_para[3];
    V_Minus(v_sat_para, v_sat, v_sat_vert, 3);

    double v_sat_para_norm = V_Norm2(v_sat_para, 3);
    double ev_sat_para[3];
    V_Divid(ev_sat_para, v_sat_para, v_sat_para_norm, 3);//

    double Vec_cross[3];
    V_Cross(Vec_cross, r_sat, ev_sat_para);
    double ev_sat_cross[3];
    V_Divid(ev_sat_cross, Vec_cross, r_norm_sat, 3);//


    const double r_target[3] = { rv_target[0], rv_target[1], rv_target[2] };
    double target_dot_sat = V_Dot(r_target, r_sat, 3);

    if (target_dot_sat < 0.0) return false;

    // 计算卫星和目标的连线到平行地面速度方向的投影
    double r_sat_target[3];
    V_Minus(r_sat_target, r_sat, r_target, 3);
    double sat_target_dot_v_sat_para = V_Dot(r_sat_target, ev_sat_para, 3);
    double mid_dist1 = fabs(sat_target_dot_v_sat_para); // �߳�һ��
    double h2 = Re_km * sqrt(1. - (mid_dist1 / Re_km) * (mid_dist1 / Re_km));
    double h_sat_center = r_norm_sat - h2;
    const double alpha1 = atan(mid_dist1 / h_sat_center);
    if (alpha1 > half_cone_angle_rad) return false;

    // 计算卫星和目标的连线到另一条曲边的投影
    double sat_target_dot_vec_cross = V_Dot(r_sat_target, ev_sat_cross, 3);
    double mid_dist2 = fabs(sat_target_dot_vec_cross);
    double h2_cross = Re_km * sqrt(1. - (mid_dist2 / Re_km) * (mid_dist2 / Re_km));
    double h_sat_center_cross = r_norm_sat - h2_cross;
    const double alpha2 = atan(mid_dist2 / h_sat_center_cross);
    if (alpha2 > half_cone_angle_rad) return false;

    return true;
}

void AccessPointObjects(
    const double rv0[6],            // 初始卫星状态（位置和速度）
    double t_start,           // 开始时间（秒）
    double t_end,             // 结束时间（秒）
    double dt,                // 时间步长（秒）
    int num_targets,          // 目标数量
    std::vector<std::vector<double>>& results // 输出：可见性结果列表
) {
    results.resize(num_targets);

    double t = t_start;
    double rv_sat[6];
    memcpy(rv_sat, rv0, 6 * sizeof(double));

    // 定义卫星的半视场角（例如 10 度，转换为弧度）
    double half_cone_angle = 19.5 * D2R; //留一些余量

    // 传播卫星到时间 t
    double rv_sat_t[6];

    memcpy(rv_sat_t, rv0, 6 * sizeof(double));

    // 循环时间步
    while (true) {
        // 对于每个目标
        for (int target_id = 0; target_id < num_targets; ++target_id)

        {
            // 如果目标已经被观测过，跳过
            // 获取目标在时间 t 的位置
            double rv_target[3];
            double Geodetic[2];
            get_target_geogetic(target_id, t, Geodetic);
            double r = V_Norm2(rv_sat_t, 3);
            double latitude_sat = asin(rv_sat_t[2] / r);

            //先判断纬度是否在范围内
            if (fabs(Geodetic[0] - latitude_sat) < 6.0 * D2R)
            {
                Geodetic2J2000(Geodetic, rv_target, t);

                //get_target_R(target_id, t, rv_target);
                // 判断可见性
                if (is_target_visible(rv_sat_t, rv_target, half_cone_angle)) {
                    // 目标可见，记录结果
                    results[target_id].push_back(t);
                }
            }
        }

        t += dt;

        if (t > t_end) break;

        // 传播卫星到时间 t
        int flag = propagate_j2(rv_sat_t, rv_sat_t, t - dt, t);
        if (flag != 1) {
            std::cerr << "Orbit propagation failed at time " << t << std::endl;
            break;
        }
    }

    return;
}


void MultiSat_AccessPointObjects(
    std::vector<std::vector<double>> rv0_list,      // 初始卫星状态（位置和速度） 多个
    double t_start,           // 开始时间（秒）
    double t_end,             // 结束时间（秒）
    double dt,                // 时间步长（秒）
    int num_targets,          // 目标数量
    std::vector<std::vector<double>>& results // 输出：可见性结果列表
) {
    results.clear();
    results.resize(num_targets);

    for (int i = 0; i < rv0_list.size(); i++)
    {
        double rv0[6];
        for (int j = 0; j < 6; j++)
        {
            rv0[j] = rv0_list[i][j];
        }

        //double* rv0 = rv0_list[i];
        std::vector<std::vector<double>> results_temp;
        AccessPointObjects(rv0, t_start, t_end, dt, num_targets, results_temp);
        for (int j = 0; j < results.size(); j++)
        {
            //把results_temp 填到results后面
            results[j].insert(results[j].end(), results_temp[j].begin(), results_temp[j].end());
        }
    }

    for (int j = 0; j < results.size(); j++)
    {
        //从小到大排序
        std::sort(results[j].begin(), results[j].end());
    }


}
