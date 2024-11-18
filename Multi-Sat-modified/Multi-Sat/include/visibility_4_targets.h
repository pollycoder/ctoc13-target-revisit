#ifndef SATELLITE_VISIBILITY_H
#define SATELLITE_VISIBILITY_H
#include <vector>
#include <cmath>
#include <iostream>
#include "J2propagation.h" // ����֮ǰ����ĺ���
#include "OrbitMath.h"
#include <iomanip>
#include <cstring>
/**
 * @brief �жϵ���Ŀ���Ƿ�����ǿɼ��� ע�⣺�ú�������Բ׶�����Ǿ��Σ������޸ģ�����
 *
 * �ú����������Ǻ�Ŀ���λ���������Լ����ǵİ�׶�ǣ��������Ŀ���Ƿ�����ǿɼ���
 * ��������������Ŀ���Ƿ��ڵ����ͬһ�ࡣ
 *
 * @param rv_sat ���ǵ�ECI����λ������������Ϊ3��6����ʹ��ǰ3��Ԫ�أ���
 * @param rv_target Ŀ���ECI����λ������������Ϊ3����
 * @param half_cone_angle_rad ���ǵİ�׶�ǣ���λ�����ȡ�
 * @return ���Ŀ������ǿɼ�������true�����򷵻�false��
 */
bool is_target_visible(const double* rv_sat, const double* rv_target, double half_cone_angle_rad);


/**
 * @brief ����������Ŀ����һ��ʱ���ڶ����ǵĿɼ��ԡ�
 *
 * �ú�������ʼʱ�� t_start ��ʼ����ʱ�䲽�� dt �����ǹ�����д�����ֱ������ʱ�� t_end��
 * ��ÿ��ʱ�䲽������ÿ������Ŀ���Ƿ�����ǿɼ�����������洢�� results �����С�
 * results ������ÿ��Ԫ�ض�Ӧһ��Ŀ�꣬������Ŀ��ɼ���ʱ���б�
 *
 * @param rv0 ���ǳ�ʼ״̬������λ�ú��ٶȣ���ECI����ϵ������Ϊ6��
 * @param t_start ��ʼʱ�䣬��λ���롣
 * @param t_end ����ʱ�䣬��λ���롣
 * @param dt ʱ�䲽������λ���롣
 * @param num_targets ����Ŀ���������
 * @param results ������������ڴ洢ÿ��Ŀ��Ŀɼ�ʱ���б�
 */
void AccessPointObjects(
    const double rv0[6],                         // ��ʼ����״̬������λ�ú��ٶȣ�
    double t_start,                              // ��ʼʱ�䣬��λ����
    double t_end,                                // ����ʱ�䣬��λ����
    double dt,                                   // ʱ�䲽������λ����
    int num_targets,                             // ����Ŀ������
    std::vector<std::vector<double>>& results,    // �����ÿ��Ŀ��Ŀɼ�ʱ���б�
    const double half_cone_angle = 20.0 * D2R
);



/**
 * @brief ������ָ������Ŀ����һ��ʱ���ڶ�ĳ�����ǵĿɼ��ԡ�
 *
 * �ú�������ʼʱ�� t_start ��ʼ����ʱ�䲽�� dt �����ǹ�����д�����ֱ������ʱ�� t_end��
 * ��ÿ��ʱ�䲽������ÿ������Ŀ���Ƿ�����ǿɼ�����������洢�� results �����С�
 * results ������ÿ��Ԫ�ض�Ӧһ��Ŀ�꣬������Ŀ��ɼ���ʱ���б�
 *
 * @param rv0 ���ǳ�ʼ״̬������λ�ú��ٶȣ���ECI����ϵ������Ϊ6��
 * @param t_start ��ʼʱ�䣬��λ���롣
 * @param t_end ����ʱ�䣬��λ���롣
 * @param dt ʱ�䲽������λ���롣
 * @param num_targets ����Ŀ���������
 * @param results ������������ڴ洢ÿ��Ŀ��Ŀɼ�ʱ���б�
 */
void AccessPointCertainObjects(
    const double rv0[6],                         // ��ʼ����״̬������λ�ú��ٶȣ�
    double t_start,                              // ��ʼʱ�䣬��λ����
    double t_end,                                // ����ʱ�䣬��λ����
    double dt,                                   // ʱ�䲽������λ����
    const std::vector<int>& target_ids,          // ָ���ĵ���Ŀ��������
    std::vector<std::vector<double>>& results,    // �����ÿ��Ŀ��Ŀɼ�ʱ���б�
    const double half_cone_angle = 20.0 * D2R
);


/**
 * @brief ���������ǶԶ������Ŀ����һ��ʱ���ڵĿɼ��ԡ�
 *
 * �ú����Զ�����ǽ��й������������ָ����ʱ����ڣ���ʱ�䲽�� dt��
 * ����ÿ�����Ƕ�ÿ������Ŀ��Ŀɼ��ԡ�����洢�� results �����У�
 * ����ÿ��Ԫ�ض�Ӧһ��Ŀ�꣬������Ŀ�걻�������ǹ۲⵽��ʱ���б�
 *
 * @param rv0_list ���ǳ�ʼ״̬������λ�ú��ٶȣ����б�����Ϊ std::vector<double*>��
 *                 ÿ��Ԫ���ǳ���Ϊ6�� double ����ָ�룬��ʾһ�����ǵĳ�ʼ״̬��
 * @param t_start ��ʼʱ�䣬��λ���롣
 * @param t_end ����ʱ�䣬��λ���롣
 * @param dt ʱ�䲽������λ���롣
 * @param num_targets ����Ŀ���������
 * @param results ������������ڴ洢ÿ��Ŀ��Ŀɼ�ʱ���б�
 *                ����Ϊ std::vector<std::vector<double>>��
 *                results[i] �洢Ŀ�� i ���������ǹ۲⵽��ʱ���б�
 */
void MultiSat_AccessPointObjects(
    std::vector<std::vector<double>> rv0_list,      // ��ʼ����״̬��λ�ú��ٶȣ� ���
    double t_start,           // ��ʼʱ�䣨�룩
    double t_end,             // ����ʱ�䣨�룩
    double dt,                // ʱ�䲽�����룩
    int num_targets,          // Ŀ������
    std::vector<std::vector<double>>& results, // ������ɼ��Խ���б�
    const double half_cone_angle = 20.0 * D2R
);


int test_visibilty();

#endif // SATELLITE_VISIBILITY_H