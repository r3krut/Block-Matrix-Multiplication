#include <vector>
#include <chrono>

#include "utils.h"

//size of matrix
#define M_SIZE 2880

int main(int argc, char *argv[])
{
    std::srand( (unsigned)time(0) );

    const size_t m_size = 1000;
    const size_t b_size = 20;

    double s_time = 0.0, in_time = 0.0, ext_time = 0.0, seq_time = 0.0;

    std::string pa = "matrices/a.txt";
    std::string pb = "matrices/b.txt";
    double **mat_a = create_matrix(m_size, 0);
    double **mat_b = create_matrix(m_size, 1);
    double **mat_c = simple_multiplication<double>(mat_a, mat_b, m_size, s_time);

    double *lin_mat_c = split_on_blocks<double>(mat_c, m_size, b_size, 2); //this is etalon for checking

    double *lin_a = split_on_blocks<double>(mat_a, m_size, b_size, 1);
    double *lin_b = split_on_blocks<double>(mat_b, m_size, b_size, 0);
    write_to_file(lin_a, m_size*(m_size+b_size)/2, pa);
    write_to_file(lin_b, m_size*(m_size+b_size)/2, pb);

//    print_lin_mtx<double>(lin_a, m_size*(m_size+b_size)/2);
//    print_lin_mtx<double>(lin_b, m_size*(m_size+b_size)/2);

//    std::cout << "\nMatrix of result(common mult. lin. repr.):\n";
//    print_lin_mtx<double>(lin_mat_c, m_size*m_size);

    std::map<std::pair<size_t, size_t>, double*> mmat_a = read_from_file<double>(m_size, b_size, pa, 1);
    std::map<std::pair<size_t, size_t>, double*> mmat_b = read_from_file<double>(m_size, b_size, pb, 0);

//    std::cout << "\n";
//    for (auto it = mmat_a.begin(); it != mmat_a.end(); it++)
//    {
//        std::pair<std::pair<size_t, size_t>, double*> m = *it;
//        std::cout << m.first.first << "|" << m.first.second << "\n";
//        print_lin_mtx<double>(m.second, b_size*b_size);
//    }

//    std::cout << "\n\n";

//    for (auto it = mmat_b.begin(); it != mmat_b.end(); it++)
//    {
//        std::pair<std::pair<size_t, size_t>, double*> m = *it;
//        std::cout << m.first.first << "|" << m.first.second << "\n";
//        print_lin_mtx<double>(m.second, b_size*b_size);
//    }

//    print_matrix<double>(mat_a, m_size);
//    print_matrix<double>(mat_b, m_size);
//    print_matrix<double>(mat_c, m_size);

    double *seq_lin_c = seq_block_mat_multiplication(mmat_a, mmat_b, m_size, b_size, seq_time);
//    std::cout << "Block mult. result:\n";
    double *ext_lin_c = external_parallel_block_mat_multiplication(mmat_a, mmat_b, m_size, b_size, 4, ext_time);
//    print_lin_mtx(ext_lin_c, m_size*m_size);
    double *int_lin_c = internal_parallel_block_mat_multiplication(mmat_a, mmat_b, m_size, b_size, 4, in_time);
//    print_lin_mtx(int_lin_c, m_size*m_size);

    release_map_matrix(mmat_a, m_size, b_size, 1);
    release_map_matrix(mmat_b, m_size, b_size, 0);

    std::cout << "\n";

    bool seq_check = compare_two_matrices(seq_lin_c, lin_mat_c, m_size);
    if (seq_check)
        std::cout << "Correct multiplication(Seq)\n";
    else
        std::cout << "Incorrect multiplication(Seq)\n";

    bool ext_check = compare_two_matrices(ext_lin_c, lin_mat_c, m_size);
    if (ext_check)
        std::cout << "Correct multiplication(External)\n";
    else
        std::cout << "Incorrect multiplication(External)\n";

    bool int_check = compare_two_matrices(int_lin_c, lin_mat_c, m_size);
    if (int_check)
        std::cout << "Correct multiplication(Internal)\n";
    else
        std::cout << "Incorrect multiplication(Internal)\n";

    std::cout << "\nTime simple mult.:         " << s_time << "\n";
    std::cout << "Time seq mult.:              " << seq_time << "\n";
    std::cout << "Time external mult.:         " << ext_time << "\n";
    std::cout << "Time internal mult.:         " << in_time << "\n";
    std::cout << "Ratio (simp./seq.):                     " << s_time / seq_time << "\n";
    std::cout << "Ratio (simp./ext.):                     " << s_time / ext_time << "\n";
    std::cout << "Ratio (simp./int.):                     " << s_time / in_time << "\n";

    return 0;
}
