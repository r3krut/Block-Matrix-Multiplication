#include <omp.h>
#include <vector>
#include <chrono>

#include "utils.h"

#include <map>

//size of matrix
#define N 10


void foo(std::string afn, std::string bfn)
{
    size_t m_size = 4;
    size_t b_size = 2;
    size_t elem_count = m_size*(m_size + b_size) / 2; //count of elems whitch store

    std::map<std::pair<int, int>, double*> A;
    double* A_elements = new double[elem_count];
    std::map<std::pair<int, int>, double*> B;
    double* B_elements = new double[elem_count];

    //read file
    std::ifstream Ainf(afn);
    std::ifstream Binf(bfn);

    double a11, a12, a21, a22;
    double b11, b12, b21, b22;
    int A_cnt = 0;
    int B_cnt = 0;

    for(int i = 1; i <= m_size / b_size; i++)
    {
        for(int j = 1; j <= i; j++)
        {
            Ainf >> a11 >> a12 >> a21 >> a22;
            A[std::pair<int, int>(i, j)] = A_elements + A_cnt;
            *(A_elements + A_cnt++) = a11;
            *(A_elements + A_cnt++) = a12;
            *(A_elements + A_cnt++) = a21;
            *(A_elements + A_cnt++) = a22;

            Binf >> b11 >> b12 >> b21 >> b22;
            B[std::pair<int, int>(j, i)] = B_elements + B_cnt;
            *(B_elements + B_cnt++) = b11;
            *(B_elements + B_cnt++) = b12;
            *(B_elements + B_cnt++) = b21;
            *(B_elements + B_cnt++) = b22;
        }
    }

    //create transparent blocks for A matrix
    for(int i = 1; i < m_size / b_size; i++)
    {
        for(int j = i + 1; j <= m_size / b_size; j++)
        {
            double* block = new double[4];
            *(block) = A[std::make_pair(j, i)][0];
            *(block + 1) = A[std::make_pair(j, i)][2];
            *(block + 2) = A[std::make_pair(j, i)][1];
            *(block + 3) = A[std::make_pair(j, i)][3];
            A[std::make_pair(i, j)] = block;
        }
    }

    double* result = new double [m_size * m_size];
    int cnt = 0;

    for(int j = 1; j <= m_size / b_size; j++)
    {
        for(int i = 1; i <= m_size / b_size; i++)
        {
            double* res_block = (result + cnt);
            cnt += b_size * b_size;
            for(int k = 1; k <= m_size / b_size; k++)
            {
                //prevent zero multiplication
                if(k > j)
                    continue;

                double* block = block_multiplication(A[std::make_pair(i , k)], B[std::make_pair(k, j)], 2);
                for(int i = 0; i < 4; i++)
                {
                    *(res_block + i) = *(res_block + i) + *(block + i);
                }
                delete [] block;
            }
        }
    }

    for(int i = 0; i < m_size * m_size; i++)
    {
        std::cout << result[i] << " ";
        if(i % 4 == 3)
            std::cout << std::endl;
    }

    Ainf.close();
    Binf.close();
}

int main(int argc, char *argv[])
{
    foo("A4.txt", "B4.txt");

    return 0;
//    std::srand( (unsigned)time(0) );

//    double **a = new double*[4];
//    double **b = new double*[4];
//    for (size_t i = 0; i < 4; i++)
//    {
//        a[i] = new double[4];
//        b[i] = new double[4];
//    }
//    fill_from_file(a, b, 4, "test.txt");

//    std::cout << "Matrix A:\n";
//    print_matrix(a, 4);
//    std::cout << "Matrix B:\n";
//    print_matrix(b, 4);
//    std::cout << "\n";

//    double *lin_a = split_on_blocks(a, 4, 2, 1); //size of block is 2. Row accomodation
//    double *lin_b = split_on_blocks(b, 4, 2, 0); //size of block is 2. Column accomodation

//    std::cout << "Lin repr. mat. A:\n";
//    print_lin_mtx<double>(lin_a, 16);
//    std::cout << "Lin repr. mat. B:\n";
//    print_lin_mtx<double>(lin_b, 16);

//    double **c = common_multiplication<double>(a, b, 4);
//    std::cout << "Matrix C:\n";
//    print_matrix(c, 4);

//    double *lin_c = block_matrix_multiplication<double>(lin_a, lin_b, 16, 2);
//    std::cout << "Lin repr. mat. C:\n";
//    print_lin_mtx<double>(lin_c, 16);

//    double **mat_a = create_matrix(N, 0);
//    double **mat_b = create_matrix(N, 0);

//    print_matrix(mat_a, N);

//    double *blin_mat_a = split_on_blocks(mat_a, N, 2, 0);
//    print_lin_mtx<double>(blin_mat_a, N * N);
//    double *bl2 = split_on_blocks(mat_a, N, 2, 0);
//    std::cout << "Matrix[1]\n";
//    print_matrix(mat_a, N);
//    std::cout << "\nMatrix[2]\n";
//    print_matrix(mat_b, N);

//    std::chrono::steady_clock::time_point c_st = std::chrono::steady_clock::now();
//    double **mat_cc = common_multiplication<double>(mat_a, mat_b, N);
//    std::chrono::steady_clock::time_point c_fn = std::chrono::steady_clock::now();
//    std::chrono::duration<double> c_time_span = std::chrono::duration_cast<std::chrono::duration<double> >(c_fn - c_st);

//    std::chrono::steady_clock::time_point p_st = std::chrono::steady_clock::now();
//    double **mat_cp = parallel_multiplication<double>(mat_a, mat_b, N, 2);
//    std::chrono::steady_clock::time_point p_fn = std::chrono::steady_clock::now();
//    std::chrono::duration<double> p_time_span = std::chrono::duration_cast<std::chrono::duration<double> >(p_fn - p_st);

//    bool check = compare_two_matrices(mat_cc, mat_cp, N);
//    if (!check)
//    {
//        std::cerr << "Matrices is not equals!\n";
//        return 0;
//    }

//    std::cout << "Time of common multiplication(Sec): " << c_time_span.count() << "\n";
//    std::cout << "Time of parallel multiplication(Sec): " << p_time_span.count() << "\n";
//    std::cout << "Rate: " << c_time_span.count() / p_time_span.count() << "\n\n";

//    release_matrix(mat_a, N);
//    release_matrix(mat_b, N);
//    release_matrix(mat_cc, N);
//    release_matrix(mat_cp, N);

    return 0;
}
