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

    double a_els[b_size * b_size];
    double b_els[b_size * b_size];
    int A_cnt = 0;
    int B_cnt = 0;

    for(int i = 1; i <= m_size / b_size; i++)
    {
        for(int j = 1; j <= i; j++)
        {
            A[std::pair<int, int>(i, j)] = A_elements + A_cnt;
            for (size_t k = 0; k < b_size * b_size; k++)
            {
                Ainf >> a_els[k];
                *(A_elements + A_cnt++) = a_els[k];
            }

            B[std::pair<int, int>(j, i)] = B_elements + B_cnt;
            for (size_t k = 0; k < b_size * b_size; k++)
            {
                Binf >> b_els[k];
                *(B_elements + B_cnt++) = b_els[k];
            }
        }
    }

    //create transpose blocks for A matrix
    for(int i = 1; i < m_size / b_size; i++)
    {
        for(int j = i + 1; j <= m_size / b_size; j++)
        {
            double* block = transpose_linear_matrix(A[std::make_pair(j, i)], b_size);
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

                double* block = block_multiplication(A[std::make_pair(i , k)], B[std::make_pair(k, j)], b_size);
                for(int l = 0; l < 4; l++)
                {
                    *(res_block + l) = *(res_block + l) + *(block + l);
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
//    foo("tests/A4.txt", "tests/B4.txt");
//    std::map<std::pair<size_t, size_t>, double*> mat_a = read_from_file<double>(4, 2, "tests/A4.txt", 1);
//    std::map<std::pair<size_t, size_t>, double*> mat_b = read_from_file<double>(4, 2, "tests/B4.txt", 0);

//    //matrices
//    for (auto it = mat_a.begin(); it != mat_a.end(); it++)
//    {
//       std::pair<std::pair<size_t, size_t>, double*> m = *it;
//       print_lin_mtx<double>(m.second, 4);
//    }

//    for (auto it = mat_b.begin(); it != mat_b.end(); it++)
//    {
//       std::pair<std::pair<size_t, size_t>, double*> m = *it;
//       print_lin_mtx<double>(m.second, 4);
//    }

//    double *mat_c = seq_block_mat_multiplication<double>(mat_a, mat_b, 4, 2);
//    print_lin_mtx<double>(mat_c, 4);

    std::srand( (unsigned)time(0) );
    const size_t m_size = 6;
    const size_t b_size = 2;
    std::string pa = "matrices/a.txt";
    std::string pb = "matrices/b.txt";
    double **mat_a = create_matrix(m_size, 0);
    double **mat_b = create_matrix(m_size, 1);
    double **mat_c = common_multiplication<double>(mat_a, mat_b, m_size);

    double *lin_mat_c = split_on_blocks<double>(mat_c, m_size, b_size, 2); //this is etalon for checking

//    print_matrix<double>(mat_a, m_size);
//    print_matrix<double>(mat_b, m_size);
//    print_matrix<double>(mat_c, m_size);
    std::cout << "\n";

    double *lin_a = split_on_blocks<double>(mat_a, m_size, b_size, 1);
    double *lin_b = split_on_blocks<double>(mat_b, m_size, b_size, 0);
    write_to_file(lin_a, m_size, b_size, pa);
    write_to_file(lin_b, m_size, b_size, pb);

    print_lin_mtx<double>(lin_a, m_size*(m_size+b_size)/2);
    print_lin_mtx<double>(lin_b, m_size*(m_size+b_size)/2);

    std::cout << "\nMatrix of result(common mult. lin. repr.):\n";
    print_lin_mtx<double>(lin_mat_c, m_size*m_size);

    std::map<std::pair<size_t, size_t>, double*> mmat_a = read_from_file<double>(m_size, b_size, pa, 1);
    std::map<std::pair<size_t, size_t>, double*> mmat_b = read_from_file<double>(m_size, b_size, pb, 0);

    std::cout << "\n";
    for (auto it = mmat_a.begin(); it != mmat_a.end(); it++)
    {
        std::pair<std::pair<size_t, size_t>, double*> m = *it;
        std::cout << m.first.first << "|" << m.first.second << "\n";
        print_lin_mtx<double>(m.second, b_size*b_size);
    }

    std::cout << "\n\n";

    for (auto it = mmat_b.begin(); it != mmat_b.end(); it++)
    {
        std::pair<std::pair<size_t, size_t>, double*> m = *it;
        std::cout << m.first.first << "|" << m.first.second << "\n";
        print_lin_mtx<double>(m.second, b_size*b_size);
    }

    print_matrix<double>(mat_a, m_size);
    print_matrix<double>(mat_b, m_size);
    print_matrix<double>(mat_c, m_size);

    std::cout << "Block mult. result:\n";
    double *lin_c = seq_block_mat_multiplication(mmat_a, mmat_b, m_size, b_size);
    print_lin_mtx(lin_c, m_size*m_size);

    return 0;

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
