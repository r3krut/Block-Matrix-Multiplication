#include "utils.h"

double dRand()
{
    return ((double)(rand() % 20000 - 10001) / 100.0);
}

float fRand()
{
    return ((float)(rand() % 20000 - 10001) / 100.0);
}

double** create_double_matrix(const size_t m_size, const bool _type)
{
    assert(m_size >= 2);

    double **mat = new double*[m_size];
    for (size_t i = 0; i < m_size; i++)
        mat[i] = new double[m_size];

    if (_type)
    {
        //top-triangular
        for (size_t i = 0; i < m_size; i++)
            for (size_t j = 0; j < m_size; j++)
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
        for (size_t i = 0; i < m_size; i++)
            for (size_t j = i; j < m_size; j++)
            {
                mat[i][j] = dRand();
                mat[j][i] = mat[i][j];
            }
    }
    return mat;
}

float** create_float_matrix(const size_t m_size, const bool _type)
{
    assert(m_size >= 2);

    float **mat = new float*[m_size];
    for (size_t i = 0; i < m_size; i++)
        mat[i] = new float[m_size];

    if (_type)
    {
        //top-triangular
        for (size_t i = 0; i < m_size; i++)
            for (size_t j = 0; j < m_size; j++)
            {
                if (i > j)
                    mat[i][j] = 0.0;
                else
                    mat[i][j] = fRand();
            }
    }
    else
    {
        //symmetric
        for (size_t i = 0; i < m_size; i++)
            for (size_t j = i; j < m_size; j++)
            {
                mat[i][j] = fRand();
                mat[j][i] = mat[i][j];
            }
    }
    return mat;
}

void create_csv_file(std::list<std::pair<size_t, double> > &list, const std::string &file_name)
{
    assert(list.size() > 0 && file_name.size() != 0);

    std::ofstream of(file_name);
    if (!of.is_open())
        throw std::runtime_error("File '" + file_name + "' don't was created.\n");

    of << "Block size" << "," << "Time(s.)\n";
    for (auto it = list.begin(); it != list.end(); it++)
        of << std::to_string((*it).first) << "," << std::to_string((*it).second) << "\n";
    of.close();
}

void make_new_dir(const std::string &dir_name)
{
    assert(dir_name.size() > 0);

    std::string command = "mkdir " + dir_name;
    size_t stat = system(command.c_str());
    if (stat)
    {
        std::string what = "Dir '" + dir_name + "' don't was created.\n";
        throw std::system_error(stat, std::system_category(), what);
    }
}
