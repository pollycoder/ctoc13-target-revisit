#include "random_threads.h"

//��[min,max)��Χ���������ʵ��
double realRand(const double& min, const double& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}
//��[min,max]��Χ�������������
int intRand(const int& min, const int& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(generator);
}