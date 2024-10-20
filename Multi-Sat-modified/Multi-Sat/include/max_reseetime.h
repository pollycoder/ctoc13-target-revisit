#include <memory.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdint>

#include "Constant.h"
#include"OrbitFun.h"
#include"OrbitMath.h"
#include <vector>
#include <algorithm>
std::vector<double> max_reseetime(std::vector<std::vector<double>> reseetime_all, double t_start, double t_end, std::vector<double>& reseetime_max);
std::vector<double> max_reseetime(std::vector<std::vector<double>> reseetime_all, std::vector<double>& reseetime_max);
std::vector<double> max_reseetime(std::vector<std::vector<uint16_t>> reseetime_all, std::vector<double>& reseetime_max);

// 已经访问过的时刻表的最大重访时间（不算最后一个gap）
// 判断是否停止扩展时可用
// 可能会出现某个目标点因为一直没被重访过导致max_gap其实为0的，这没有关系
// 最终每个节点都会被看到，一定会有一个非空的时刻表，只要这个时间超过了，也可以停止扩展
void max_revisit_interval_beforeEnd(std::vector<double>& reseetime_max,const std::vector<std::vector<double>>& reseetime_all);