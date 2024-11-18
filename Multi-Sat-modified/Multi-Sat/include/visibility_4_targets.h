#ifndef SATELLITE_VISIBILITY_H
#define SATELLITE_VISIBILITY_H
#include <vector>
#include <cmath>
#include <iostream>
#include "J2propagation.h" // 包含之前定义的函数
#include "OrbitMath.h"
#include <iomanip>
#include <cstring>
/**
 * @brief 判断地面目标是否对卫星可见。 注意：该函数仅用圆锥，不是矩形，可以修改！！！
 *
 * 该函数根据卫星和目标的位置向量，以及卫星的半锥角，计算地面目标是否对卫星可见。
 * 它考虑了卫星与目标是否在地球的同一侧。
 *
 * @param rv_sat 卫星的ECI坐标位置向量（长度为3或6，仅使用前3个元素）。
 * @param rv_target 目标的ECI坐标位置向量（长度为3）。
 * @param half_cone_angle_rad 卫星的半锥角，单位：弧度。
 * @return 如果目标对卫星可见，返回true；否则返回false。
 */
bool is_target_visible(const double* rv_sat, const double* rv_target, double half_cone_angle_rad);


/**
 * @brief 计算多个地面目标在一段时间内对卫星的可见性。
 *
 * 该函数从起始时间 t_start 开始，以时间步长 dt 对卫星轨道进行传播，直到结束时间 t_end。
 * 在每个时间步，计算每个地面目标是否对卫星可见，并将结果存储在 results 向量中。
 * results 向量的每个元素对应一个目标，包含该目标可见的时间列表。
 *
 * @param rv0 卫星初始状态向量（位置和速度），ECI坐标系，长度为6。
 * @param t_start 起始时间，单位：秒。
 * @param t_end 结束时间，单位：秒。
 * @param dt 时间步长，单位：秒。
 * @param num_targets 地面目标的数量。
 * @param results 输出参数，用于存储每个目标的可见时间列表。
 */
void AccessPointObjects(
    const double rv0[6],                         // 初始卫星状态向量（位置和速度）
    double t_start,                              // 起始时间，单位：秒
    double t_end,                                // 结束时间，单位：秒
    double dt,                                   // 时间步长，单位：秒
    int num_targets,                             // 地面目标数量
    std::vector<std::vector<double>>& results,    // 输出：每个目标的可见时间列表
    const double half_cone_angle = 20.0 * D2R
);



/**
 * @brief 计算多个指定地面目标在一段时间内对某颗卫星的可见性。
 *
 * 该函数从起始时间 t_start 开始，以时间步长 dt 对卫星轨道进行传播，直到结束时间 t_end。
 * 在每个时间步，计算每个地面目标是否对卫星可见，并将结果存储在 results 向量中。
 * results 向量的每个元素对应一个目标，包含该目标可见的时间列表。
 *
 * @param rv0 卫星初始状态向量（位置和速度），ECI坐标系，长度为6。
 * @param t_start 起始时间，单位：秒。
 * @param t_end 结束时间，单位：秒。
 * @param dt 时间步长，单位：秒。
 * @param num_targets 地面目标的数量。
 * @param results 输出参数，用于存储每个目标的可见时间列表。
 */
void AccessPointCertainObjects(
    const double rv0[6],                         // 初始卫星状态向量（位置和速度）
    double t_start,                              // 起始时间，单位：秒
    double t_end,                                // 结束时间，单位：秒
    double dt,                                   // 时间步长，单位：秒
    const std::vector<int>& target_ids,          // 指定的地面目标编号向量
    std::vector<std::vector<double>>& results,    // 输出：每个目标的可见时间列表
    const double half_cone_angle = 20.0 * D2R
);


/**
 * @brief 计算多个卫星对多个地面目标在一段时间内的可见性。
 *
 * 该函数对多个卫星进行轨道传播，并在指定的时间段内，以时间步长 dt，
 * 计算每个卫星对每个地面目标的可见性。结果存储在 results 向量中，
 * 其中每个元素对应一个目标，包含该目标被任意卫星观测到的时间列表。
 *
 * @param rv0_list 卫星初始状态向量（位置和速度）的列表，类型为 std::vector<double*>，
 *                 每个元素是长度为6的 double 数组指针，表示一个卫星的初始状态。
 * @param t_start 起始时间，单位：秒。
 * @param t_end 结束时间，单位：秒。
 * @param dt 时间步长，单位：秒。
 * @param num_targets 地面目标的数量。
 * @param results 输出参数，用于存储每个目标的可见时间列表。
 *                类型为 std::vector<std::vector<double>>，
 *                results[i] 存储目标 i 被任意卫星观测到的时间列表。
 */
void MultiSat_AccessPointObjects(
    std::vector<std::vector<double>> rv0_list,      // 初始卫星状态（位置和速度） 多个
    double t_start,           // 开始时间（秒）
    double t_end,             // 结束时间（秒）
    double dt,                // 时间步长（秒）
    int num_targets,          // 目标数量
    std::vector<std::vector<double>>& results, // 输出：可见性结果列表
    const double half_cone_angle = 20.0 * D2R
);


int test_visibilty();

#endif // SATELLITE_VISIBILITY_H