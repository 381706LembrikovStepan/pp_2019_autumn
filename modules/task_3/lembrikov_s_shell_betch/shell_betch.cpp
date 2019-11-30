// Copyright 2019 Lembrikov Stepan
#include <../../../modules/task_3/lembrikov_s_shell_betch/shell_betch.h>
#include <algorithm>
#include <vector>
#include <iomanip>

std::vector<int> getOtrVector(int n) {
    std::vector <int> a(n);
    int k = 1;
    for (int i = 0; i < n; i++) {
        k = -k;
        a[i] = 1 * k;
    }
    return a;
}

std::vector<int> getRandomVector(int n) {
    // std::mt19937 engine;
    std::vector<int> a(n);
    std::mt19937 gen(time(0));
    // std::uniform_real_distribution<> urd(0, 1);
    std::uniform_int_distribution <int> urd(0, n);
    // engine.seed(n);
    // int k = 1;
    for (int i = 0; i < n; i++) {
        // k = -k;
        a[i] = urd(gen);
    }
    return a;
}

std::vector <int> Chet_Betch(const std::vector <int> &mas_1,
    const std::vector <int> &mas_2) {
    int size1 = mas_1.size();
    int size2 = mas_2.size();
    int size_res = size1 / 2 + size2 / 2 + size1 % 2 + size2 % 2;
    std::vector <int> mas_res(size_res);
    int it1 = 0;
    int it2 = 0;
    int i = 0;

    while ((it1 < size1) && (it2 < size2)) {
        if (mas_1[it1] <= mas_2[it2]) {
            mas_res[i] = mas_1[it1];
            it1 += 2;
        } else {
            mas_res[i] = mas_2[it2];
            it2 += 2;
        }
        i++;
    }

    if (it1 >= size1) {
        for (int j = it2; j < size2; j += 2) {
            mas_res[i] = mas_2[j];
            i++;
        }
    } else {
        for (int j = it1; j < size1; j += 2) {
            mas_res[i] = mas_1[j];
            i++;
        }
    }
    return mas_res;
}

std::vector <int> Nechet_Betch(const std::vector <int> &mas_1,
    const std::vector <int> &mas_2) {
    int size1 = mas_1.size();
    int size2 = mas_2.size();
    int size_res = size1 / 2 + size2 / 2;
    std::vector <int> mas_res(size_res);
    int it1 = 1;
    int it2 = 1;
    int i = 0;

    while ((it1 < size1) && (it2 < size2)) {
        if (mas_1[it1] <= mas_2[it2]) {
            mas_res[i] = mas_1[it1];
            it1 += 2;
        } else {
            mas_res[i] = mas_2[it2];
            it2 += 2;
        }
        i++;
    }

    if (it1 >= size1) {
        for (int j = it2; j < size2; j += 2) {
            mas_res[i] = mas_2[j];
            i++;
        }
    } else {
        for (int j = it1; j < size1; j += 2) {
            mas_res[i] = mas_1[j];
            i++;
        }
    }
    return mas_res;
}

std::vector <int> Sravnenie_Chet_Nechet(const std::vector <int> &mas_res_1,
    const std::vector <int> &mas_res_2) {
    int size1 = mas_res_1.size();
    int size2 = mas_res_2.size();
    int flag = 0;
    if (size1 - size2 == 2)
        flag = 1;
    if (size2 - size1 == 2)
        flag = 2;
    int size_res = mas_res_1.size() + mas_res_2.size();
    int size_min = size_res / 2 - 1;
    std::vector <int> mas_res(size_res);
    int buf;
    int i = 0;
    if (flag == 0) {
        for (i = 0; i < size1; i++) {
            mas_res[2 * i] = mas_res_1[i];
            mas_res[2 * i + 1] = mas_res_2[i];
        }
    }

    if ((flag == 1) || (flag == 2)) {
        for (i = 0; i < size_min; i++) {
            mas_res[2 * i] = mas_res_1[i];
            mas_res[2 * i + 1] = mas_res_2[i];
        }
    }

    if (flag == 1) {
        mas_res[2 * i] = mas_res_1[i];
        mas_res[2 * i + 1] = mas_res_1[i + 1];
    }

    if (flag == 2) {
        mas_res[2 * i] = mas_res_2[i];
        mas_res[2 * i + 1] = mas_res_2[i + 1];
    }

    for (int i = 1; i < size_res; i++) {
        if (mas_res[i] < mas_res[i - 1]) {
            buf = mas_res[i - 1];
            mas_res[i - 1] = mas_res[i];
            mas_res[i] = buf;
        }
    }

    return mas_res;
}

std::vector <int> ShellSort(const std::vector <int> &mas, int size) {
    int step, i, j, tmp;
    std::vector <int> mas_res(mas);
    for (step = size / 2; step > 0; step /= 2)
        for (i = step; i < size; i++)
            for (j = i - step; j >= 0 && mas_res[j] > mas_res[j + step]; j -= step) {
                tmp = mas_res[j];
                mas_res[j] = mas_res[j + step];
                mas_res[j + step] = tmp;
            }
    return mas_res;
}


std::vector <int> Shell(std::vector <int> mas) {
    int size;
    int rank;
    int ost;
    // int k;
    int flag = 0;
    int ostatok = 0;
    int k = 0;
    int size_mas = mas.size();
    if (size_mas == 1)
        return mas;
    int it_proizved = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ((size % 2 == 1) && (size > 1)) {
        flag = -1;
        size = size - 1;
    } else if ((size % 2 == 1) && (size == 1)) {
        flag = 1;
    }
    ostatok = size_mas % size;
    k = size_mas / size;
    std::vector <int> part_mas(k + ostatok, 0);
    MPI_Status status;

    for (int i = 0; i < size - 2; i += 2) {
        if (rank == i) {
            part_mas = ShellSort({ mas.cbegin() + i * k, mas.cbegin() + (i + 1) * k }, k);
            std::copy(part_mas.begin(), part_mas.begin() + k, mas.begin() + i * k);
            if (i < size - 1) {
                MPI_Sendrecv(&part_mas[0], k, MPI_INT, i + 1, 0,
                    &mas[(i + 1) * k], k, MPI_INT, i + 1, 0, MPI_COMM_WORLD, &status);
            }
        }
    }
    for (int i = 1; i < size - 1; i += 2) {
        if (rank == i) {
            part_mas = ShellSort({ mas.cbegin() + i * k, mas.cbegin() + (i + 1) * k }, k);
            std::copy(part_mas.begin(), part_mas.begin() + k, mas.begin() + k * i);
            MPI_Sendrecv(&part_mas[0], k, MPI_INT, i - 1, 0,
                &mas[(i - 1) * k], k, MPI_INT, i - 1, 0, MPI_COMM_WORLD, &status);
        }
    }

    if (rank == size - 2) {
        part_mas = ShellSort({ mas.cend() - 2 * k - ostatok, mas.cend() - k - ostatok}, k);
        std::copy(part_mas.begin(), part_mas.begin() + k, mas.end() - 2 * k - ostatok);
        MPI_Sendrecv(&part_mas[0], k, MPI_INT, size - 1, 0,
            &mas[(size - 1) * k], k + ostatok, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &status);
    }

    if (rank == size - 1) {
        part_mas = ShellSort({ mas.cend() - k - ostatok, mas.cend() }, k + ostatok);
        std::copy(part_mas.begin(), part_mas.begin() + k + ostatok, mas.end() - k - ostatok);
        MPI_Sendrecv(&part_mas[0], k + ostatok, MPI_INT, size - 2, 0,
            &mas[(size - 2) * k], k, MPI_INT, size - 2, 0, MPI_COMM_WORLD, &status);
    }

    int flag_nechet_proc = 0;
    int iter = 0;
    int count_iter = 0;
    int buf = size;
    while (buf != 1) {
        if (buf % 2 == 0) {
            buf = buf / 2;
            count_iter++;
        } else {
            buf--;
            count_iter++;
        }
    }
    if (size > 1) {
        while (iter < count_iter) {
        if ((k % 2) == 1)
            ost = 1;
        else
            ost = 0;
        std::vector <int> res_part_mas1(k + ost);
        std::vector <int> res_part_mas2(k - ost);
        std::vector <int> res_part_mas(2 * k);
        int it_par = 0;
        int smesh = 0;
        int dobavka = 0;
        int it = (size_mas - ostatok) / k / 2;
        if (flag_nechet_proc == 1) {
            it++;
        }

        for (int i = 0; i < it; i++) {
            smesh = it_par * 2 * k;

            if (i == (it - 1)) {
                if (flag_nechet_proc == 1) {
                    ost = (k / 2) % 2;
                    res_part_mas1.resize(k / 2 + k / 4 + ost + ostatok / 2 + ostatok % 2);
                    res_part_mas2.resize(k / 2 + k / 4 + ostatok / 2);
                    res_part_mas.resize(k + k / 2 + ostatok);
                    dobavka = k / 2 - k + ostatok;
                }
                if ((flag_nechet_proc == 0) && ((size_mas - ostatok) == (k * 2 * it))) {
                    if (ost == 0) {
                        res_part_mas1.resize(k + ostatok / 2 + ostatok % 2);
                        res_part_mas2.resize(k + ostatok / 2);
                    }
                    if (ost == 1) {
                        res_part_mas1.resize(k + ost + ostatok / 2);
                        res_part_mas2.resize(k - ost + ostatok / 2 + ostatok % 2);
                    }
                    res_part_mas.resize(2 * k + ostatok);
                    dobavka = ostatok;
                }
            }
            if (flag_nechet_proc == 1) {
                flag_nechet_proc = 0;
            }

            if (rank == it_par * 2) {
                res_part_mas1 = Chet_Betch({ mas.cbegin() + smesh,
    mas.cbegin() + smesh + k },
                    { mas.cbegin() + smesh + k, mas.cbegin() + smesh + 2 * k  + dobavka });
                MPI_Send(&res_part_mas1[0], res_part_mas1.size(), MPI_INT, it_par * 2 + 1, 0, MPI_COMM_WORLD);
            }
            if (rank == it_par * 2 + 1) {
                MPI_Status status0;
                res_part_mas2 = Nechet_Betch({ mas.cbegin() + smesh,
    mas.cbegin() + smesh + k },
                    { mas.cbegin() + smesh + k, mas.cbegin() + smesh + 2 * k + dobavka });
                MPI_Recv(&res_part_mas1[0], res_part_mas1.size(), MPI_INT, it_par * 2, 0, MPI_COMM_WORLD, &status0);
                res_part_mas = Sravnenie_Chet_Nechet(res_part_mas1, res_part_mas2);
                if (flag == 0) {
                    MPI_Send(&res_part_mas[0], res_part_mas.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
                } else if (flag == -1) {
                    MPI_Send(&res_part_mas[0], res_part_mas.size(), MPI_INT, size, 1, MPI_COMM_WORLD);
                }
            }
            if (flag == 0) {
                if (rank == 0) {
                    MPI_Status status1;
                    if (size > 1) {
                        MPI_Recv(&mas[smesh], res_part_mas.size(),
                            MPI_INT, it_par * 2 + 1, 1, MPI_COMM_WORLD, &status1);
                    }
                }
                MPI_Bcast(&mas[smesh], res_part_mas.size(), MPI_INT, 0, MPI_COMM_WORLD);
            } else if (flag == -1) {
                if (rank == size) {
                    MPI_Status status1;
                    if (size > 1) {
                        MPI_Recv(&mas[smesh], res_part_mas.size(),
                            MPI_INT, it_par * 2 + 1, 1, MPI_COMM_WORLD, &status1);
                    }
                }
                MPI_Bcast(&mas[smesh], res_part_mas.size(), MPI_INT, size, MPI_COMM_WORLD);
            }
            it_par++;
        }
        if ((((size_mas - ostatok) / k) % 2) != 0) {
            flag_nechet_proc = 1;
        }
        k = k * 2;
        iter++;
    }
    }
    return mas;
}
