#include "utils.h"

double dRand()
{
    return ((double)(rand() % 10 ));
}

double** create_matrix(const size_t n, const bool _type)
{
    assert(n >= 2);

    double **mat = new double*[n];
    for (size_t i = 0; i < n; i++)
        mat[i] = new double[n];

    if (_type)
    {
        //top-triangular
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
            {
                if (i > j)
                    mat[i][j] = 0.0;
                else
                    mat[i][j] = dRand();
            }
    }
    else
    {
        //symmetric
        for (size_t i = 0; i < n; i++)
            for (size_t j = i; j < n; j++)
            {
                mat[i][j] = dRand();
                mat[j][i] = mat[i][j];
            }
    }
    return mat;
}

void release_matrix(double **mat, const size_t n)
{
    for (size_t i = 0; i < n; i++)
        delete [] mat[i];
    delete [] mat;

    std::cout << "Matrix was released.\n";
}

void print_matrix(double **mat, const size_t n)
{
    assert(n >= 2);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            std::cout << mat[i][j] << " ";
        std::cout << "\n";
    }
    std::cout << "Size of matrix: " << n << "x" << n << "\n";
}

void write_to_file(const double *lin_mat, const size_t n, const std::string &file_name)
{
    std::ofstream of(file_name);
    if (!of.is_open())
        throw std::runtime_error("File '" + file_name + "' dot't was created.\n");
    for (size_t i = 0; i < n; i++)
        of << lin_mat[i] << " ";
    of.close();
}

double *mtx_to_linear(double **mat, const size_t n)
{
    assert(n >= 2 && mat != NULL);

    double *lin_mtx = new double[n * n];
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            lin_mtx[i*n+j] = mat[i][j];
    return lin_mtx;
}

double* split_on_blocks(double **mat, const size_t n, const size_t block_sz, const bool _type)
{
    assert(mat != NULL && n >= 2 && block_sz <= n && n % block_sz == 0);

    double *lin_repr = new double[n * n];
    size_t lin_ind = 0;
    size_t count = 0;
    if (_type) //by rows
    {
        for (size_t i = 0; i < n; i += block_sz)
            for (size_t j = 0; j < n; j += block_sz)
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
    else //by columns
    {
        for (size_t j = 0; j < n; j += block_sz)
            for (size_t i = 0; i < n; i += block_sz)
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

    return lin_repr;
}

void fill_from_file(double **a, double **b, const size_t n, const std::string &fn)
{
    std::ifstream inf(fn);

    double aa, bb;
    size_t i = 0, j = 0;
    while (inf >> aa >> bb)
    {
        a[i][j] = aa;
        b[i][j] = bb;
        j++;
        if (j == 4)
        {
            j = 0;
            i += 1;
        }
    }
    inf.close();
}
