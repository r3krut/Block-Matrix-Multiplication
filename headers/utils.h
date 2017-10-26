#ifndef UTILS_H
#define UTILS_H

#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <exception>
#include <assert.h>
#include <omp.h>

//generate double numbers function
double dRand();

/**
 * @brief create_matrix - perform creation of matrix 'mat'
 * @param n
 * @param _type - this is type of matrix: 0 - symmetric, 1 - top-triangular
 * @return double pointer
 */
double** create_matrix(const size_t n, const bool _type);
void release_matrix(double **mat, const size_t n);

/**
 * @brief split_on_blocks - performs splitting on blocks and returns the matrix in the linearal representation
 * @param mat
 * @param n
 * @param block_sz - size of block
 * @param _type - type of accomodation. If 1 - accomodation by rows, by columns
 * @return
 */
double* split_on_blocks(double **mat, const size_t n, const size_t block_sz, const bool _type);

void write_to_file(const double *lin_mat, const size_t n, const std::string &file_name);

template<typename ELEM_T>
void read_from_file(ELEM_T *lin_mat, const size_t size, const std::string file_name)
{
    std::ifstream inf(file_name);
    if (!inf.is_open())
        throw std::runtime_error("File '" + file_name + "' dot't was created.\n");

    ELEM_T elem;
    size_t ind = 0;
    while (inf >> elem)
        lin_mat[ind++] = static_cast<ELEM_T>(elem);
    inf.close();
}

/**
 * @brief mtx_to_linear - transform source matrix to linearal representation
 * @param mat
 * @param n
 * @return pointer on a doubles
 */
double* mtx_to_linear(double **mat, const size_t n);

template<typename LMTX_T>
void release_linear_mtx(LMTX_T *lin_mtx)
{
    delete [] lin_mtx;
}

template<typename MTX_T>
bool compare_two_matrices(MTX_T **mat_a, MTX_T **mat_b, const size_t n)
{
    assert(mat_a != NULL && mat_b != NULL && n >= 2);

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (mat_a[i][j] != mat_b[i][j])
                return false;
    return true;
}

template<typename MTX_T>
/**
 * @brief common_multiplication - multiplication of two matrices
 * @param mat_a
 * @param mat_b
 * @param n
 * @return new matrix
 */
MTX_T** common_multiplication(MTX_T **mat_a, MTX_T **mat_b, const size_t n)
{
    assert(mat_a != NULL && mat_b != NULL && n >= 2);

    MTX_T **mat_c = new MTX_T*[n];
    for (size_t i = 0; i < n; i++)
        mat_c[i] = new MTX_T[n];

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
        {
            mat_c[i][j] = 0.0;
            for (size_t k = 0; k < n; k++)
                mat_c[i][j] += mat_a[i][k] * mat_b[k][j];
        }

    return mat_c;
}

template<typename MTX_T>
/**
 * @brief parallel_multiplication - simple parallel multiplication without splitting on blocks
 * @param mat_a
 * @param mat_b
 * @param n
 * @param num_thrds
 * @return
 */
MTX_T** parallel_multiplication(MTX_T **mat_a, MTX_T **mat_b, const size_t n, const int num_thrds)
{
    assert(mat_a != NULL && mat_b != NULL && n >= 2);

    MTX_T **mat_c = new MTX_T*[n];
    for (size_t i = 0; i < n; i++)
        mat_c[i] = new MTX_T[n];

#pragma omp parallel for num_threads(num_thrds)
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
        {
            mat_c[i][j] = 0.0;
            for (size_t k = 0; k < n; k++)
                mat_c[i][j] += mat_a[i][k] * mat_b[k][j];
        }

    return mat_c;
}

template<typename BLK_T>
/**
 * @brief block_multiplication - multiplication of two linearal block
 *
 *       |3 1 4 3|         |5 3 2 0|
 *       |2 4 5 3|         |4 6 0 1|
 *   A = |1 5 3 6|    B =  |6 6 2 3|
 *       |7 7 4 3|         |3 4 7 9|
 *
 *                                  For block size 2
 *
 * @param block1. Example: 3124-first block, 4353-second block, etc. This is the block-row represenation, i.e. 31244353...
 * @param block2. Example: 5346-first block, 6634-second block, etc. This is the block-column representation, i.e. 53466634..
 * @param block_sz
 * @return new block - result of multiplication. It's also a block representation
 */
BLK_T* block_multiplication(BLK_T *block1, BLK_T *block2, const size_t block_sz)
{
    assert(block_sz >= 1);

    size_t len_of_block = block_sz * block_sz;
    BLK_T *res_block = new BLK_T[len_of_block];
    size_t ind_res_blk = 0;
    for (size_t i = 0; i < len_of_block; i += block_sz) //traverse through all the rows in the first block. Step is 'block_sz'
    {
        for (size_t j = 0; j < block_sz; j++)
        {
            size_t k = j;
            BLK_T sum = 0.0;
            for (size_t ii = i; ii < ((i + block_sz) >= len_of_block ? len_of_block : (i + block_sz)); ii++)
            {
                sum += block1[ii] * block2[k];
                k += block_sz;
            }
            res_block[ind_res_blk++] = sum;
        }
    }
    return res_block;
}

template<typename MTX_T>
MTX_T* sum_two_blocks(MTX_T *block1, MTX_T *block2, const size_t block_sz)
{
    assert(block_sz >= 1);

    MTX_T *res_block = new MTX_T[block_sz * block_sz];
    for (size_t i = 0; i < block_sz * block_sz; i++)
        res_block[i] = block1[i] + block2[i];
    return res_block;
}

template<typename MTX_T>
/**
 * @brief seq_block_mat_multiplication - performs multiplication of two linear matrices which represented as a block matrices
 * @param lmat_a
 * @param lmat_b
 * @param n - size of matrix. Size of linear matrix is n * n
 * @param block_sz - size of block
 * @return
 */
MTX_T* seq_block_mat_multiplication(MTX_T *lmat_a, MTX_T *lmat_b, const size_t n, const size_t block_sz)
{
    assert(lmat_a != NULL && lmat_b != NULL && n % block_sz == 0 && n >= 2 && block_sz >= 1);

    size_t elem_inblock = block_sz * block_sz;
    size_t lin_mat_size = n * n;
    size_t count_blocks = lin_mat_size / elem_inblock;
    size_t count_rows = n / block_sz; //count of blocks equals the count of rows

    MTX_T *lmat_c = new MTX_T[n * n];

    return lmat_c;
}

template<typename LMTX_T>
void print_lin_mtx(const LMTX_T *lin_mtx, const size_t n)
{
    assert (lin_mtx != NULL && n >= 2);

    for (size_t i = 0; i < n; i++)
        std::cout << lin_mtx[i] << " ";
    std::cout << "\nSize of linear matrix: " << n << "\n";
}
void print_matrix(double **mat, const size_t n);

void fill_from_file(double **a, double **b, const size_t n, const std::string &fn);

#endif // UTILS_H
