// Copyright 2019 Lembrikov Stepan

#include <mpi.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>
#include <../../../modules/task_1/lembrikov_s_min_elem_vector/min_elem_vector.h>

std::vector<int> getIdentityVector(int n) {
    std::vector <int> a(n);
    for (int i = 0; i < n; i++) {
        a[i] = 1;
    }
    return a;
}

std::vector<int> getRandomVector(int n) {
    std::mt19937 engine;
    engine.seed(n);
    std::vector<int> a(n);
    for (int i = 0; i < n; i++) {
        if (engine() % 100 >= 0) {
            a[i] = engine() % 100;
        } else {
            a[i] = 0;
        }
    }
    return a;
}

std::vector<int> getConstVector(int n, int c) {
    std::vector <int> a(n);
    for (int i = 0; i < n; i++) {
        a[i] = c;
    }
    return a;
}

int MinOfVector(const std::vector <int> a, int n) {
    int res = 0;
    int size;
    int rank;
//    std::vector <int> i = getRandomVector(7);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            if (i % size) {
                MPI_Send(&a[i], 1, MPI_INT, i%size, 0, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Status status;
    int b;
    int prom_res;
    if (rank == 0) {
        for (int i = 0; i < n; i += size) {
            prom_res = std::min(a[i], prom_res);
        }
    } else {
            for (int i = rank; i < n; i += size) {
                MPI_Recv(&b, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
                prom_res = std::min(b, prom_res);
            }
        }
    MPI_Reduce(&prom_res, &res, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    return res;
}