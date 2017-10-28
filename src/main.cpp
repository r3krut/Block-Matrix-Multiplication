#include <vector>
#include <chrono>

#include "utils.h"

//size of matrix
#define N 10

int main(int argc, char *argv[])
{
    std::srand( (unsigned)time(0) );

    const size_t m_size = 4;
    const size_t b_size = 2;
    std::string pa = "matrices/a.txt";
    std::string pb = "matrices/b.txt";
    double **mat_a = create_matrix(m_size, 0);
    double **mat_b = create_matrix(m_size, 1);
    double **mat_c = common_multiplication<double>(mat_a, mat_b, m_size);

    double *lin_mat_c = split_on_blocks<double>(mat_c, m_size, b_size, 2); //this is etalon for checking

    double *lin_a = split_on_blocks<double>(mat_a, m_size, b_size, 1);
    double *lin_b = split_on_blocks<double>(mat_b, m_size, b_size, 0);
    write_to_file(lin_a, m_size*(m_size+b_size)/2, pa);
    write_to_file(lin_b, m_size*(m_size+b_size)/2, pb);

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

    release_map_matrix(mmat_a, m_size, b_size, 1);
    release_map_matrix(mmat_b, m_size, b_size, 0);

    bool check = compare_two_matrices(lin_c, lin_mat_c, m_size);
    if (check)
        std::cout << "Correct multiplication\n";
    else
        std::cout << "Incorrect multiplication\n";

    std::cout << "\n\n";

    auto it = mmat_b.begin();
    for (size_t i = 0; i < m_size*(m_size+b_size)/2; i++)
    {
        std::cout << (*it).second[i] << " ";
    }
    std::cout << "\n";

    it = mmat_a.begin();
    for (size_t i = 0; i < m_size*m_size; i++)
    {
        std::cout << (*it).second[i] << " ";
    }
    std::cout << "\n";

    return 0;
}
