#include "visibility_4_targets.h"
#include "OrbitFun.h"

#include <algorithm>

// 计算两个向量之间的距离
//double distance(const double& x1, const double& y1, const double& z1,
//    const double& x2, const double& y2, const double& z2) {
//    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
//}
//
//// 计算两个向量之间的方位角（绕Z轴旋转的角度，单位：弧度）
//double azimuth(const double& x, const double& y) {
//    return atan2(y, x);
//}
//
//// 计算两个向量之间的仰角（绕X轴旋转的角度，单位：弧度）
//double elevation(const double& x, const double& y, const double& z) {
//    return atan2(z, sqrt(x * x + y * y));
//}

// 判断点是否在四棱锥视场内
//bool is_target_visible(const double* rv_sat, const double* r_target, double half_cone_angle_rad) {
//    // 提取探测点和被观测点的位置速度
//    double rv_sat_x = rv_sat[0];
//    double rv_sat_y = rv_sat[1];
//    double rv_sat_z = rv_sat[2];
//    double r_target_x = r_target[0];
//    double r_target_y = r_target[1];
//    double r_target_z = r_target[2];
//
//    // 计算目标点到地心的距离
//    double distance_to_target = distance(0.0, 0.0, 0.0, r_target_x, r_target_y, r_target_z);
//
//    // 计算观测点的方位角和仰角
//    double azimuth_sat = azimuth(rv_sat_x, rv_sat_y);
//    double elevation_sat = elevation(rv_sat_x, rv_sat_y, rv_sat_z);
//
//    // 计算目标点相对于观测点的方向向量
//    double direction_vector_x = r_target_x - rv_sat_x;
//    double direction_vector_y = r_target_y - rv_sat_y;
//    double direction_vector_z = r_target_z - rv_sat_z;
//
//    // 计算目标点相对于观测点的方位角和仰角
//    double target_azimuth = azimuth(direction_vector_x, direction_vector_y);
//    double target_elevation = elevation(direction_vector_x, direction_vector_y, direction_vector_z);
//
//    // 定义四棱锥的半锥角（单位：弧度）
//    double half_cone_angle = half_cone_angle_rad;
//
//    // 判断目标点是否在四棱锥视场内
//    if (abs(target_azimuth - azimuth_sat) <= half_cone_angle && target_elevation >= elevation_sat - half_cone_angle) {
//        return true;
//    }
//
//    return false;
//}

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


// 矩形视角的可见角度平方
double angle_target_visible(const double* rv_sat, const double* rv_target, double half_cone_angle_rad) {
    const double r_sat[3] = { rv_sat[0], rv_sat[1], rv_sat[2] };
    const double v_sat[3] = { rv_sat[3], rv_sat[4], rv_sat[5] };
    const double r_dot_v_sat = V_Dot(r_sat, v_sat, 3);
    const double r_norm_sat = V_Norm2(r_sat, 3);
    const double v_norm_sat = V_Norm2(v_sat, 3);

    if (r_norm_sat <= Re_km) // 卫星位置不合理
        return 10000.0;

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

    if (target_dot_sat < 0.0) return 10000.0;

    double result_angle = 0.0;

    // 计算卫星和目标的连线到平行地面速度方向的投影
    double r_sat_target[3];
    V_Minus(r_sat_target, r_sat, r_target, 3);
    double sat_target_dot_v_sat_para = V_Dot(r_sat_target, ev_sat_para, 3);
    double mid_dist1 = fabs(sat_target_dot_v_sat_para); // �߳�һ��
    if (mid_dist1 >= Re_km) {
        // 说明目标水平位移太大，直接返回不可见 or 返回一个极大值
        return 10000.0;
    }
    double h2 = Re_km * sqrt(1. - (mid_dist1 / Re_km) * (mid_dist1 / Re_km));
    double h_sat_center = r_norm_sat - h2;
    const double alpha1 = atan(mid_dist1 / h_sat_center);
    if (alpha1 > half_cone_angle_rad) result_angle += (alpha1 - half_cone_angle_rad)* (alpha1 - half_cone_angle_rad);

    // 计算卫星和目标的连线到另一条曲边的投影
    double sat_target_dot_vec_cross = V_Dot(r_sat_target, ev_sat_cross, 3);
    double mid_dist2 = fabs(sat_target_dot_vec_cross);
    if (mid_dist2 >= Re_km) {
        // 说明目标水平位移太大，直接返回不可见 or 返回一个极大值
        return 10000.0;
    }
    double h2_cross = Re_km * sqrt(1. - (mid_dist2 / Re_km) * (mid_dist2 / Re_km));
    double h_sat_center_cross = r_norm_sat - h2_cross;
    const double alpha2 = atan(mid_dist2 / h_sat_center_cross);
    if (alpha2 > half_cone_angle_rad) result_angle += (alpha2 - half_cone_angle_rad) * (alpha2 - half_cone_angle_rad);

    return result_angle;
}


double angle_earth_view_angle(const double* rv_sat, const double* rv_target)
{
    double dot = V_Dot(rv_sat, rv_target, 3);
    return acos(dot / V_Norm2(rv_sat, 3) / V_Norm2(rv_target, 3));
}


void AccessPointObjects(
    const double rv0[6],            // 初始卫星状态（位置和速度）
    double t_start,           // 开始时间（秒）
    double t_end,             // 结束时间（秒）
    double dt,                // 时间步长（秒）
    int num_targets,          // 目标数量
    std::vector<std::vector<double>>& results, // 输出：可见性结果列表
    const double half_cone_angle
) {
    results.resize(num_targets);

    double t;
    
    double rv_sat[6];
    memcpy(rv_sat, rv0, 6 * sizeof(double));

    if (fmod(t_start, dt) < 1.0e-10) {
        t = t_start;
    }
    else {
        // t不是60的倍数，先进它的下一个60s
        t = t_start - fmod(t_start, dt) + dt;
        propagate_j2(rv_sat, rv_sat, t_start, t);
        //int flag;
        //rv02rvf(flag, rv_sat, rv_sat, t - t_start, mu_km_s);
    }

    // 定义卫星的半视场角（例如 10 度，转换为弧度）
    //double half_cone_angle = 20.0 * D2R; //留一些余量

    // 传播卫星到时间 t
    double rv_sat_t[6];

    memcpy(rv_sat_t, rv_sat, 6 * sizeof(double));
    

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
            //if (fabs(Geodetic[0] - latitude_sat) < 6.0 * D2R)
            //{
                Geodetic2J2000(Geodetic, rv_target, t);

                //get_target_R(target_id, t, rv_target);
                // 判断可见性
                if (is_target_visible(rv_sat_t, rv_target, half_cone_angle)) {
                    // 目标可见，记录结
                        results[target_id].push_back(t);
                }
           // }
        }

        t += dt;

        if (t > t_end) break;
        if (t > 2.0 * 86400.0) break;

        // 传播卫星到时间 t
        int flag = propagate_j2(rv_sat_t, rv_sat_t, t - dt, t);
        //rv02rvf(flag, rv_sat_t, rv_sat_t, dt, mu_km_s);
        if (flag != 1) {
            //std::cerr << "Orbit propagation failed at time " << t << std::endl;
            break;
        }
    }

    return;
}

void AccessPointObjects_linearJ2 (
    const double rv0[6],            // 初始卫星状态（位置和速度）
    double t_start,           // 开始时间（秒）
    double t_end,             // 结束时间（秒）
    double dt,                // 时间步长（秒）
    int num_targets,          // 目标数量
    std::vector<std::vector<double>>& results, // 输出：可见性结果列表
    const double half_cone_angle
) {
    results.resize(num_targets);

    double t;

    double rv_sat[6];
    memcpy(rv_sat, rv0, 6 * sizeof(double));

    if (fmod(t_start, dt) < 1.0e-10) {
        t = t_start;
    }
    else {
        // t不是60的倍数，先进它的下一个60s
        t = t_start - fmod(t_start, dt) + dt;
        propagate_linearJ2(rv_sat, rv_sat, t_start, t);
        //int flag;
        //rv02rvf(flag, rv_sat, rv_sat, t - t_start, mu_km_s);
    }

    // 定义卫星的半视场角（例如 10 度，转换为弧度）
    //double half_cone_angle = 20.0 * D2R; //留一些余量

    // 传播卫星到时间 t
    double rv_sat_t[6];

    memcpy(rv_sat_t, rv_sat, 6 * sizeof(double));


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
            //if (fabs(Geodetic[0] - latitude_sat) < 6.0 * D2R)
            //{
            Geodetic2J2000(Geodetic, rv_target, t);

            //get_target_R(target_id, t, rv_target);
            // 判断可见性
            if (is_target_visible(rv_sat_t, rv_target, half_cone_angle)) {
                // 目标可见，记录结
                results[target_id].push_back(t);
            }
            // }
        }

        t += dt;

        if (t > t_end) break;
        if (t > 2.0 * 86400.0) break;

        // 传播卫星到时间 t
        propagate_linearJ2(rv_sat_t, rv_sat_t, t - dt, t);
        //rv02rvf(flag, rv_sat_t, rv_sat_t, dt, mu_km_s);
        //if (flag != 1) {
        //    //std::cerr << "Orbit propagation failed at time " << t << std::endl;
        //    break;
        //}
    }

    return;
}



void AccessPointCertainObjects(
    const double rv0[6],                         // 初始卫星状态向量（位置和速度）
    double t_start,                              // 起始时间，单位：秒
    double t_end,                                // 结束时间，单位：秒
    double dt,                                   // 时间步长，单位：秒
    const int& target_id,          // 指定的地面目标编号向量
    std::vector<double>& results,    // 输出：每个目标的可见时间列表
    const double half_cone_angle
) {
    results.resize(0);

    double t = t_start;
    double rv_sat[6];
    memcpy(rv_sat, rv0, 6 * sizeof(double));

    // 定义卫星的半视场角（例如 10 度，转换为弧度）
    //double half_cone_angle = 20.0 * D2R; //留一些余量

    // 传播卫星到时间 t
    double rv_sat_t[6];

    memcpy(rv_sat_t, rv0, 6 * sizeof(double));

    // 循环时间步
    while (true) {
        // 对于每个目标

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
                results.push_back(t);
            }
        }


        t += dt;

        if (t > t_end) break;
        if (t > 2.0 * 86400.0) break;

        // 传播卫星到时间 t
        int flag = propagate_j2(rv_sat_t, rv_sat_t, t - dt, t);
        //rv02rvf(flag, rv_sat_t, rv_sat_t, dt, mu_km_s);
        if (flag != 1) {
            //std::cerr << "Orbit propagation failed at time " << t << std::endl;
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
    std::vector<std::vector<double>>& results, // 输出：可见性结果列表
    const double half_cone_angle
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
        AccessPointObjects(rv0, t_start, t_end, dt, num_targets, results_temp, half_cone_angle);
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

//简化版：只用60s步长
void AccessPeriodCertainObject(
    const double* rv0,
    const double& t_standard,
    const double* dv,
    const int& target_id,
    double* t_period,
    const double& half_cone_angle
) {
    const double tf_m = t_standard - fmod(t_standard, 60.0);
    const double tf_p = tf_m + 60.0;
    double rv_temp_m[6], rv_temp_p[6], rv0_copy[6];
    double step_m = 60.0, step_p = 60.0;

    memcpy(rv0_copy, rv0, 6 * sizeof(double));
    propagate_j2(rv0_copy, rv_temp_m, t_standard, tf_m);

    for (int i = 0; i < 3; i++) {
        rv0_copy[i + 3] += dv[i];
    }
    propagate_j2(rv0_copy, rv_temp_p, t_standard, tf_p);

    //先查两个检测点
    bool ifVisible_m, ifVisible_p;

    double target_r_m[3], target_r_p[3];
    get_target_R(target_id, tf_m, target_r_m);
    ifVisible_m = is_target_visible(rv_temp_m, target_r_m, half_cone_angle);
    get_target_R(target_id, tf_p, target_r_p);
    ifVisible_p = is_target_visible(rv_temp_p, target_r_p, half_cone_angle);

    // ATK无法检测，放弃
    if (!ifVisible_m && !ifVisible_p) {
        t_period[0] = 1.0e10;
        t_period[1] = 1.0e10;
        return;
    }
    
    if (!ifVisible_m) t_period[0] = tf_p;
    if (!ifVisible_p) t_period[1] = tf_m;
    
    //如果头可见，从头开始向前递推
    if (ifVisible_m) {
        t_period[0] = tf_m;                                     // 从head开始向前，给定初始时刻
        double tm_temp = tf_m;
        while (ifVisible_m) {
            propagate_j2(rv_temp_m, rv_temp_m, tm_temp, tm_temp - 60.0);
            get_target_R(target_id, tm_temp - 60.0, target_r_m);
            ifVisible_m = is_target_visible(rv_temp_m, target_r_m, half_cone_angle);
            if (ifVisible_m) {
                tm_temp -= 60.0;
                if(tm_temp >= 0.0) t_period[0] = tm_temp;
            }
        }
    }
   

    //如果尾可见，从尾开始向后递推
    if (ifVisible_p) {
        t_period[1] = tf_p;                                     // 从tail开始向后，给定末端时刻
        double tp_temp = tf_p;
        while (ifVisible_p) {
            propagate_j2(rv_temp_p, rv_temp_p, tp_temp, tp_temp + 60.0);
            get_target_R(target_id, tp_temp + 60.0, target_r_p);
            ifVisible_p = is_target_visible(rv_temp_p, target_r_p, half_cone_angle);
            if (ifVisible_p) {
                tp_temp += 60.0;
                t_period[1] = tp_temp;
            }
        }
    }

}
