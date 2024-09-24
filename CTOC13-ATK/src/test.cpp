#include <iostream>

#include "ResultData.h"



int main() {
	ResultData* result = new ResultData();
	result->read_data("result.txt");
	result->write_atk();

	return 0;
}