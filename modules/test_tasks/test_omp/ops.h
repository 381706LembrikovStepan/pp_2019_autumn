// Copyright 2018 Nesterov Alexander
#ifndef MODULES_TEST_TASKS_OMP_OPS_H_
#define MODULES_TEST_TASKS_OMP_OPS_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(size_t sz);
int getParallelOperations(std::vector<int> vec, std::string ops);
int getSequentialOperations(std::vector<int> vec, std::string ops);

#endif  // MODULES_TEST_TASKS_OMP_OPS_H_