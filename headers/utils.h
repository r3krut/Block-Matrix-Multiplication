#ifndef UTILS_H
#define UTILS_H

/*                        Matrix utils file by /rekrut/
 *
 *
 * Contains functions to work with matrices. Simple multiplication,
 * several type of parallization of matrix multiplication with a help OpenMP technology and other.
 * Also this file contains functions such as read, write matrices to the file and from file and print it.
 *
 *
 * */

#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <exception>
#include <list>
#include <map>
#include <cmath>
#include <vector>
#include <omp.h>

//generate double numbers function
double dRand();
//generate float numbers
float fRand();

/**
 * @brief create_double_matrix - perform creation of matrix
 * @param n
 * @param _type - this is type of matrix: 0 - symmetric, 1 - top-triangular
 * @return double pointer
 */
double** create_double_matrix(const size_t m_size, const bool _type);
float** create_float_matrix(const size_t m_size, const bool _type);

template<typename MTX_T>
void release_matrix(MTX_T **&mat, const size_t m_size)
{
    for (size_t i = 0; i < m_size; i++)
        delete [] mat[i];
    delete [] mat;

    std::cout << "Matrix was released.\n";
}

template<typename LMTX_T>
void release_linear_mtx(LMTX_T *&lin_mtx)
{
    delete [] lin_mtx;
}

template<typename MTX_T>
/**
 * @brief release_map_matrix - release a linear matrix throught pointers which contains in the map
 * @param m
 * @param m_size
 * @param b_size
 * @param _type - type of matix. 1 - symmetric, 0 - top-triangular
 */
void release_map_matrix(std::map<std::pair<size_t, size_t>, MTX_T*> &m,
                        const size_t m_size,
                        const size_t b_size,
                        bool _type)
{
    assert(m.size() > 0 && m_size >= 2 && b_size >= 1 && m_size % b_size == 0);

    if (!_type)
    {
        auto it = m.begin();
        delete [] (*it).second;
    }
    else
    {
        for (size_t i = 1; i < m_size / b_size; i++)
        {
            //delete all transpose blocks
            for(size_t j = i + 1; j <= m_size / b_size; j++)
            {
                delete [] m[std::make_pair(i, j)];
            }
        }
        //delete blocks which placed below a main diagonal
        auto it = m.begin();
        delete [] (*it).second;
    }
}

template<typename MTX_T>
/**
 * @brief split_on_blocks - performs splitting on blocks and returns the matrix in the linearal representation
 * @param mat
 * @param n
 * @param block_sz - size of block
 * @param _type - type of accomodation. If 1 - accomodation of matrix by rows(Matrix must be symmetrical)
 *                                      If 0 - accomodation of matrix by columns(Matrix must be top-triangular)
 *                                      If 2 - accomodation of matrix by columns. It is common matrix. Needed to checking
 * @return
 */
MTX_T* split_on_blocks(MTX_T **mat, const size_t n, const size_t block_sz, const size_t _type)
{
    assert(mat != NULL && n >= 2 && block_sz <= n && n % block_sz == 0);

    size_t count_elems = n * (n + block_sz) / 2;
    MTX_T *lin_repr = new MTX_T[count_elems];
    size_t lin_ind = 0;
    size_t count = 0;
    if (_type == 1) //by rows
    {
        for (size_t i = 0; i < n; i += block_sz)
            for (size_t j = 0; j <= i; j += block_sz)
            {
                size_t inn = (i + block_sz) >= n ? n : (i + block_sz);
                size_t jnn = (j + block_sz) >= n ? n : (j + block_sz);
                for (size_t ii = i; ii < inn; ii++)
                {
                    for (size_t jj = j; jj < jnn; jj++)
                        lin_repr[lin_ind++] = mat[ii][jj];
                }
                count++;
            }
        std::cout << "Count of blocks: " << count << "\n";
    }
    else if (_type == 0) //by columns
    {
        for (size_t j = 0; j < n; j += block_sz)
            for (size_t i = 0; i <= j; i += block_sz)
            {
                size_t inn = (i + block_sz) >= n ? n : (i + block_sz);
                size_t jnn = (j + block_sz) >= n ? n : (j + block_sz);
                for (size_t ii = i; ii < inn; ii++)
                {
                    for (size_t jj = j; jj < jnn; jj++)
                        lin_repr[lin_ind++] = mat[ii][jj];
                }
                count++;
            }
        std::cout << "Count of blocks: " << count << "\n";
    }
    else if (_type == 2) //when we want convert commont matrix to a linear representation
    {
        MTX_T *lin_repr_comm_mtx = new MTX_T[n * n]; //allocate new memmory
        delete [] lin_repr; // delete previously allocated memmory
        for (size_t j = 0; j < n; j += block_sz)
            for (size_t i = 0; i < n; i += block_sz)
            {
                size_t inn = (i + block_sz) >= n ? n : (i + block_sz);
                size_t jnn = (j + block_sz) >= n ? n : (j + block_sz);
                for (size_t ii = i; ii < inn; ii++)
                {
                    for (size_t jj = j; jj < jnn; jj++)
                        lin_repr_comm_mtx[lin_ind++] = mat[ii][jj];
                }
                count++;
            }
        std::cout << "Count of blocks: " << count << "\n";
        return lin_repr_comm_mtx;
    }
    return lin_repr;
}

template<typename MTX_T>
/**
 * @brief transpose_linear_matrix - calculates a tranpose matrix by a linear matrix
 * @param lin_block - matrix shown as a linear representation
 * @param block_sz - size of block. It is size of the quadric matrix
 * @return - tranpose linear matrix
 */
inline MTX_T* transpose_linear_matrix(MTX_T *lin_block, const size_t block_sz)
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

template<typename ELEM_T>
/**
 * @brief read_from_file - read linear matrix from file
 * @param m_size - size of matrix
 * @param b_size - size of block
 * @param file_name
 * @param _type - 1 if matrix is symmetric, 0 - if matrix is a top-triangular
 * @return - map. Where: Key is a pair of indexes of a block of matrix.
 *                       Value it is a pointer to the beginning of the corresponding block in the linear representation
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
            if (_type)
                mat[std::make_pair(i, j)]= mat_elems + m_count; //for symmetric matrix
            else
                mat[std::make_pair(j, i)] = mat_elems + m_count; //for top-triangular

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
/**
 * @brief compare_two_matrices - compares two linear matrices.
 * @param lmat_a
 * @param lmat_b
 * @param m_size - size of matices. lmat_a and lmat_b must have a same size
 * @return
 */
inline bool compare_two_matrices(MTX_T *lmat_a, MTX_T *lmat_b, const size_t m_size, const double &eps)
{
    assert(lmat_a != NULL && lmat_b != NULL && m_size >= 2);

    for (size_t i = 0; i < m_size; i++)
        if (std::fabs(lmat_a[i] - lmat_b[i]) > eps)
            return false;
    return true;
}

template<typename MTX_T>
/**
 * @brief simple_multiplication - this is common multiplication of two matrices
 * @param mat_a
 * @param mat_b
 * @param n
 * @param time - this is reference variable which will be contains the time of calculation
 * @return new matrix
 */
MTX_T** simple_multiplication(MTX_T **mat_a, MTX_T **mat_b, const size_t n, double &time)
{
    assert(mat_a != NULL && mat_b != NULL && n >= 2);

    MTX_T **mat_c = new MTX_T*[n];
    for (size_t i = 0; i < n; i++)
        mat_c[i] = new MTX_T[n];

    std::chrono::steady_clock::time_point st = std::chrono::steady_clock::now();
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
        {
            mat_c[i][j] = 0.0;
            for (size_t k = 0; k < n; k++)
                mat_c[i][j] += mat_a[i][k] * mat_b[k][j];
        }
    std::chrono::steady_clock::time_point fn = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(fn - st).count();

    return mat_c;
}

template<typename MTX_T>
/**
 * @brief parallel_multiplication - simple parallel multiplication without splitting on blocks
 * @param mat_a
 * @param mat_b
 * @param n
 * @param num_thrds
 * @param time
 * @return
 */
MTX_T** parallel_multiplication(MTX_T **mat_a, MTX_T **mat_b, const size_t n, const int num_thrds, double &time)
{
    assert(mat_a != NULL && mat_b != NULL && n >= 2);

    MTX_T **mat_c = new MTX_T*[n];
    for (size_t i = 0; i < n; i++)
        mat_c[i] = new MTX_T[n];

    std::chrono::steady_clock::time_point st = std::chrono::steady_clock::now();
#pragma omp parallel for num_threads(num_thrds)
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
        {
            mat_c[i][j] = 0.0;
            for (size_t k = 0; k < n; k++)
                mat_c[i][j] += mat_a[i][k] * mat_b[k][j];
        }
    std::chrono::steady_clock::time_point fn = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(fn - st).count();

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
 * @param enable_omp - bool variable. If 1 then the function performs in the parallel mode. 0 - sequential mode
 * @param num_thrds - number of threads
 * @return new block - result of multiplication. It's also a block representation
 */
inline BLK_T* block_multiplication(BLK_T *block1, BLK_T *block2, const size_t block_sz, bool enable_omp, const size_t num_thrds)
{
    assert(block_sz >= 1);

    size_t len_of_block = block_sz * block_sz;
    BLK_T *res_block = new BLK_T[len_of_block];
#pragma omp parallel if(enable_omp) num_threads(num_thrds)
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
            res_block[i + j] = sum;
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
 * @param time
 * @return - The result returns as a linear representation of matrix where blocks placed by columns
 */
ELEM_T* seq_block_mat_multiplication(std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_a,
                                     std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_b,
                                     const size_t m_size, const size_t b_size, double &time)
{
    assert(mat_a.size() != 0 && mat_b.size() != 0 && m_size >= 2 && b_size >= 1 && m_size % b_size == 0);

    ELEM_T *mat_c = new ELEM_T[m_size * m_size];
    for (size_t i = 0; i < m_size * m_size; i++)
        mat_c[i] = 0.0;

    size_t c_count = 0;

    std::chrono::steady_clock::time_point st = std::chrono::steady_clock::now();
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

                ELEM_T *block = block_multiplication<ELEM_T>(mat_a[std::make_pair(i , k)], mat_b[std::make_pair(k, j)], b_size, 0, 1);
                for (size_t l = 0; l < b_size * b_size; l++)
                    *(res_block + l) = *(res_block + l) + *(block + l);
                delete [] block;
            }
        }
    }
    std::chrono::steady_clock::time_point fn = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(fn - st).count();

    return mat_c;
}

template<typename ELEM_T>
/**
 * @brief inner_parallel_block_mat_multiplication - this is parallel block matrix multiplication.
 *                                                  parallelization of inner cycle. Multiplying of two blocks
 * @param mat_a
 * @param mat_b
 * @param m_size
 * @param b_size
 * @param num_thrds
 * @param time
 * @return - The result returns as a linear representation of matrix where blocks placed by columns
 */
ELEM_T* internal_parallel_block_mat_multiplication(std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_a,
                                                   std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_b,
                                                   const size_t m_size, const size_t b_size, const size_t num_thrds, double &time)
{
    assert(mat_a.size() != 0 && mat_b.size() != 0 && m_size >= 2 && b_size >= 1 && m_size % b_size == 0);

    ELEM_T *mat_c = new ELEM_T[m_size * m_size];
    for (size_t i = 0; i < m_size * m_size; i++)
        mat_c[i] = 0.0;

    size_t c_count = 0;

    std::chrono::steady_clock::time_point st = std::chrono::steady_clock::now();
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

                ELEM_T *block = block_multiplication<ELEM_T>(mat_a[std::make_pair(i , k)], mat_b[std::make_pair(k, j)], b_size, 1, num_thrds);
                for (size_t l = 0; l < b_size * b_size; l++)
                    *(res_block + l) = *(res_block + l) + *(block + l);
                delete [] block;
            }
        }
    }
    std::chrono::steady_clock::time_point fn = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(fn - st).count();

    return mat_c;
}

template<typename ELEM_T>
/**
 * @brief external_parallel_block_mat_multiplication - this is parallel block matrix multiplication.
 *                                                     parallelization of external cycle.
 *                                                     Two different blocks multiply on the different kernels
 * @param mat_a
 * @param mat_b
 * @param m_size
 * @param b_size
 * @param num_thrds
 * @param time
 * @return - The result returns as a linear representation of matrix where blocks placed by columns
 */
ELEM_T* external_parallel_block_mat_multiplication(std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_a,
                                                   std::map<std::pair<size_t, size_t>, ELEM_T*> &mat_b,
                                                   const size_t m_size, const size_t b_size, const size_t num_thrds,
                                                   double &time)
{
    assert(mat_a.size() != 0 && mat_b.size() != 0 && m_size >= 2 && b_size >= 1 && m_size % b_size == 0);

    ELEM_T *mat_c = new ELEM_T[m_size * m_size];
    for (size_t i = 0; i < m_size * m_size; i++)
        mat_c[i] = 0.0;

    size_t c_count = 0;

    std::chrono::steady_clock::time_point st = std::chrono::steady_clock::now();
    for (size_t j = 1; j <= m_size / b_size; j++)
    {
        for (size_t i = 1; i <= m_size / b_size; i++)
        {
            ELEM_T *res_block = (mat_c + c_count);
            c_count += b_size * b_size;
#pragma omp parallel for num_threads(num_thrds) //two different blocks multiply on a different kernels
            for (size_t k = 1; k <= m_size / b_size; k++)
            {
                //prevent zero multiplication
                if (k > j)
                   continue;

                ELEM_T *block = block_multiplication<ELEM_T>(mat_a[std::make_pair(i , k)], mat_b[std::make_pair(k, j)], b_size, 0, 1);
                for (size_t l = 0; l < b_size * b_size; l++)
                    *(res_block + l) = *(res_block + l) + *(block + l);
                delete [] block;
            }
        }
    }
    std::chrono::steady_clock::time_point fn = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(fn - st).count();

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
/**
 * @brief print_lin_mtx - prints the matrix which shown as linear representation
 * @param lin_mtx
 * @param m_size - matrix size - it is source size
 * @param b_size - size of block
 */
void print_lin_mtx(const LMTX_T *lin_mtx, const size_t m_size)
{
    assert (lin_mtx != NULL);

    for (size_t i = 0; i < m_size; i++)
        std::cout << lin_mtx[i] << " ";
    std::cout << "\nSize of linear matrix: " << m_size << "\n";
}

template<typename MTX_T>
void print_matrix(MTX_T **mat, const size_t m_size)
{
    assert(mat != NULL && m_size >= 2);

    for (size_t i = 0; i < m_size; i++)
    {
        for (size_t j = 0; j < m_size; j++)
            std::cout << mat[i][j] << " ";
        std::cout << "\n";
    }
    std::cout << "Size of matrix: " << m_size << "x" << m_size << "\n";
}

template<typename LMTX_T>
/**
 * @brief write_to_file - write linear matrix to file
 * @param lin_mat - linear matrix
 * @param elems_count - real number of elements in the linear representation
 * @param file_name
 */
void write_to_file(const LMTX_T *lin_mat, const size_t elems_count, const std::string &file_name)
{
    assert(lin_mat != NULL && elems_count > 0);

    std::ofstream of(file_name);
    if (!of.is_open())
        throw std::runtime_error("File '" + file_name + "' dot't was created.\n");
    for (size_t i = 0; i < elems_count; i++)
        of << lin_mat[i] << " ";
    of.close();
}

template<typename LMTX_T>
LMTX_T* read_etalon_from_file(const size_t m_size, const std::string &file_name)
{
    assert(m_size >= 2 && file_name.size() > 0);

    LMTX_T *l_mat = new LMTX_T[m_size * m_size];
    LMTX_T elem;

    std::ifstream inf(file_name);
    if (!inf.is_open())
        throw std::runtime_error("File '" + file_name + "' don't was openned.");
    for (size_t i = 0; i < m_size * m_size; i++)
    {
        inf >> elem;
        l_mat[i] = static_cast<LMTX_T>(elem);
    }
    inf.close();
    return l_mat;
}

/**
 * @brief create_csv_file - create a csv file which contains follow
 * @param list
 * @param file_name
 */
void create_csv_file(std::list<std::pair<size_t, double>> &list, const std::string &file_name);
void make_new_dir(const std::string &dir_name);

#endif // UTILS_H
