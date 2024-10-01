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