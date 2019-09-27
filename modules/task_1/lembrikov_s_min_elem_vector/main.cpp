// Copyright 2019 Lembrikov Stepan

#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "../../../modules/task_1/lembrikov_s_min_elem_vector/min_elem_vector.h"

TEST(Vector_Min_MPI, Test_On_Vector_EQ) {
    std::vector <int> a(10);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a[0] = 2;
    a[1] = 1;
    a[2] = 3;
    a[3] = 0;
    a[0] = 9;
    a[1] = 8;
    a[2] = 7;
    a[3] = 6;
    a[4] = 5;
    a[5] = 4;
    a[6] = 3;
    a[7] = 2;
    a[8] = 1;
    a[9] = 0;
    int min = MinOfVector(a, 10);
    if (rank == 0) {
        EXPECT_EQ(0, min);
    }
}

TEST(Vector_Min_MPI, Test_On_Vector_NE) {
    std::vector <int> a(4);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a[0] = 2;
    a[1] = 1;
    a[2] = 3;
    a[3] = 5;
    int min = MinOfVector(a, 4);
    if (rank == 0) {
        EXPECT_NE(0, min);
    }
}

TEST(Vector_Min_MPI, test_on_identity_vector_EQ) {
    std::vector <int> a = getIdentityVector(4);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a[0] = 0;
    int min = MinOfVector(a, 4);
    if (rank == 0) {
        EXPECT_EQ(0, min);
    }
}

TEST(Vector_Min_MPI, test_on_identity_vector_NE) {
    std::vector <int> a = getIdentityVector(10);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int min = MinOfVector(a, 10);
    if (rank == 0) {
        EXPECT_NE(0, min);
    }
}

TEST(Vector_Min_MPI, test_on_random_vector_EQ) {
    std::vector <int> a = getRandomVector(4);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a[0] = 0;
    int min = MinOfVector(a, 4);
    if (rank == 0) {
        EXPECT_EQ(0, min);
    }
}

TEST(Vector_Min_MPI, test_on_random_vector_NE) {
    std::vector <int> a = getRandomVector(10);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int min = MinOfVector(a, 10);
    if (rank == 0) {
        EXPECT_NE(-1, min);
    }
}

TEST(Vector_Min_MPI, test_on_const_vector_EQ) {
    std::vector <int> a = getConstVector(5, 10);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a[0] = 0;
    int min = MinOfVector(a, 5);
    if (rank == 0) {
        EXPECT_EQ(0, min);
    }
}

TEST(Vector_Min_MPI, test_on_const_vector_NE) {
    std::vector <int> a = getConstVector(5, 10);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int min = MinOfVector(a, 5);
    if (rank == 0) {
        EXPECT_NE(0, min);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}