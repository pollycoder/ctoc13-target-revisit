#include "random_threads.h"

//在[min,max)范围内生成随机实数
double realRand(const double& min, const double& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}
//在[min,max]范围内生成随机整数
int intRand(const int& min, const int& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(generator);
}