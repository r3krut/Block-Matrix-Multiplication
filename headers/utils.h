#ifndef UTILS_H
#define UTILS_H

#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <map>
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

/**
 * @brief mtx_to_linear - transform source matrix to linearal representation
 * @param mat
 * @param n
 * @return pointer on a doubles
 */
double* mtx_to_linear(double **mat, const size_t n);
void release_matrix(double **mat, const size_t n);

template<typename LMTX_T>
void release_linear_mtx(LMTX_T *lin_mtx)
{
    delete [] lin_mtx;
}

/**
 * @brief split_on_blocks - performs splitting on blocks and returns the matrix in the linearal representation
 * @param mat
 * @param n
 * @param block_sz - size of block
 * @param _type - type of accomodation. If 1 - accomodation by rows, by columns
 * @return
 */
double* split_on_blocks(double **mat, const size_t n, const size_t block_sz, const bool _type);

template<typename MTX_T>
/**
 * @brief transpose_linear_matrix - calculates a tranpose matrix by a linear matrix
 * @param lin_block - matrix whown as a linear representation
 * @param block_sz
 * @return tranpose linear matrix
 */
MTX_T* transpose_linear_matrix(MTX_T *lin_block, const size_t block_sz)
{
    assert(lin_block != NULL && block_sz >= 1);

    MTX_T *t_lin_bock = new MTX_T[block_sz * block_sz];

    for (size_t i = 0; i < block_sz; i++)
        for (size_t j = 0; j < block_sz; j++)
        {
            if (i != j)
                t_lin_bock[j*block_sz + i] = lin_block[i*block_sz + j];
            else
                t_lin_bock[i*block_sz + i] = lin_block[i*block_sz + i];
        }
    return t_lin_bock;
}

void write_to_file(const double *lin_mat, const size_t n, const std::string &file_name);

template<typename ELEM_T>
/**
 * @brief read_from_file - read matrix from file
 * @param m_size - size of matrix
 * @param b_size - size of block
 * @param file_name
 * @param _type - 1 if matrix is symmetric, 0 - if matrix is a top-triangular
 * @return - map where: Key is a pair of indexes of a block of matrix.
 *                      Value it is a pointer to the beginning of the corresponding block in the linear representation
 */
std::map<std::pair<size_t, size_t>, ELEM_T*> read_from_file(const size_t m_size, const size_t b_size, const std::string &file_name, bool _type)
{
    assert(m_size >= 2 && b_size >= 1
           && m_size % b_size == 0
           && file_name.size() > 0);

    size_t elems_count = m_size*(m_size + b_size) / 2; //count of elements which will be stored
    std::map<std::pair<size_t, size_t>, ELEM_T*> mat;
    ELEM_T* mat_elems = new ELEM_T[elems_count];

    std::ifstream inf(file_name);
    if (!inf.is_open())
        throw std::runtime_error("File '" + file_name + "' dot't was opened.\n");

    //linear representaion of matrix
    ELEM_T m_els[b_size * b_size];
    size_t m_count = 0;

    for (size_t i = 1; i <= m_size / b_size; i++)
    {
        for (size_t j = 1; j <= i; j++)
        {
            //saving of a pointer to the beginning of block
            mat[std::make_pair(i, j)] = mat_elems + m_count;

            //reading by blocks. Size of block is b_size, then number of elems in linear
            //representation equals b_size * b_size
            for (size_t k = 0; k < b_size * b_size; k++)
            {
                inf >> m_els[k];
                *(mat_elems + m_count++) = static_cast<ELEM_T>(m_els[k]);
            }
        }
    }

    if (_type) //for symmetric matrix
    {
        //create transpose blocks for A matrix
        for (size_t i = 1; i < m_size / b_size; i++)
        {
            for(size_t j = i + 1; j <= m_size / b_size; j++)
            {
                ELEM_T* block = transpose_linear_matrix<ELEM_T>(mat[std::make_pair(j, i)], b_size);
                mat[std::make_pair(i, j)] = block;
            }
        }
    }
    inf.close();
    return mat;
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

template<typename ELEM_T>
/**
 * @brief seq_block_mat_multiplication - this is sequential method of block matrix multiplication
 * @param mat_a
 * @param mat_b
 * @param m_size
 * @param b_size
 * @return - The result returns as a linear representation of matrix where blocks placed by columns
 */
ELEM_T* seq_block_mat_multiplication(std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_a,
                                     std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_b,
                                     const size_t m_size, const size_t b_size)
{
    assert(mat_a.size() != 0 && mat_b.size() != 0 && m_size >= 2 && b_size >= 1 && m_size % b_size == 0);

    ELEM_T *mat_c = new ELEM_T[m_size * m_size];
    size_t c_count = 0;

    for (size_t j = 1; j <= m_size / b_size; j++)
    {
        for (size_t i = 1; i <= m_size / b_size; i++)
        {
            ELEM_T *res_block = (mat_c + c_count);
            c_count += b_size * b_size;
            for (size_t k = 1; k <= m_size / b_size; k++)
            {
                //prevent zero multiplication
                if (k > j)
                    continue;
                ELEM_T *block = block_multiplication<ELEM_T>(mat_a[std::make_pair(i , k)], mat_b[std::make_pair(k, j)], b_size);
                for (size_t l = 0; l < b_size * b_size; l++)
                {
                    *(res_block + l) = *(res_block + l) + *(block + l);
                }
                delete [] block;
            }
        }
    }

    return mat_c;
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

template<typename LMTX_T>
void print_lin_mtx(const LMTX_T *lin_mtx, const size_t n)
{
    assert (lin_mtx != NULL && n >= 2);

    for (size_t i = 0; i < n; i++)
        std::cout << lin_mtx[i] << " ";
    std::cout << "\nSize of linear matrix: " << n << "\n";
}
void print_matrix(double **mat, const size_t n);

#endif // UTILS_H
