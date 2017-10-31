#include <vector>
#include <chrono>
#include <iomanip>
#include <iostream>

#include "utils.h"

//possible sizes for block
const std::vector<size_t> sizes = {1, 6, 10, 15, 20, 24, 30, 36, 40, 60, 72, 80,
                                   96, 120, 144, 160, 180, 240, 360, 480, 720};
static constexpr double eps = 2.0;

/**
 * @brief generate_files - generates a linear matrices from the source matrix and writes it to the file for each block
 *                         Also generates the etalon linear matrix for compare from source matrix for each block
 */
void generate_files(const size_t &m_size)
{
    double time = 0;
    double **sym_mat = create_double_matrix(m_size, 0);
    double **toptriang_mat = create_double_matrix(m_size, 1);
    double **mult_m = parallel_multiplication<double>(sym_mat, toptriang_mat, m_size, 4, time);

    for (size_t i = 0; i < sizes.size(); i++)
    {
        std::string path = "matrices/" + std::to_string(sizes[i]);
        make_new_dir(path);

        double *lin_sym_mat = split_on_blocks<double>(sym_mat, m_size, sizes[i], 1);
        double *lin_toptriang_mat = split_on_blocks<double>(toptriang_mat, m_size, sizes[i], 0);
        double *lin_mult_mat = split_on_blocks<double>(mult_m, m_size, sizes[i], 2);

        write_to_file<double>(lin_sym_mat, m_size*(m_size+sizes[i])/2, path + "/" + "a.txt");
        write_to_file<double>(lin_toptriang_mat, m_size*(m_size+sizes[i])/2, path + "/" + "b.txt");
        write_to_file<double>(lin_mult_mat, m_size*m_size, path + "/" + "etalon.txt");

        release_linear_mtx<double>(lin_sym_mat);
        release_linear_mtx<double>(lin_toptriang_mat);
        release_linear_mtx<double>(lin_mult_mat);
    }

    release_matrix<double>(sym_mat, m_size);
    release_matrix<double>(toptriang_mat, m_size);
    release_matrix<double>(mult_m, m_size);
}

/**
 * @brief perform_float_test - test for float type
 * @param m_size
 * @param num_th
 * @return
 */
std::tuple<std::list<std::pair<size_t, double> >,
           std::list<std::pair<size_t, double> >,
           std::list<std::pair<size_t, double> > > perform_float_test(const size_t &m_size, const size_t num_th = 4)
{
    assert(m_size >= 2);

    std::list<std::pair<size_t, double> > res_seq, res_int, res_ext;
    for (size_t i = 1; i < sizes.size(); i++)
    {
        double time_seq = 0.0, time_internal = 0.0, time_external = 0.0;
        std::string path_a = "/home/rekrut/QTProjects/MultiplicationBlockMatrix/matrices/" + std::to_string(sizes[i]) + "/a.txt";
        std::string path_b = "/home/rekrut/QTProjects/MultiplicationBlockMatrix/matrices/" + std::to_string(sizes[i]) + "/b.txt";
        std::string path_etalon = "/home/rekrut/QTProjects/MultiplicationBlockMatrix/matrices/" + std::to_string(sizes[i]) + "/etalon.txt";
        std::map<std::pair<size_t, size_t>, float*> mmat_a = read_from_file<float>(m_size, sizes[i], path_a, 1);
        std::map<std::pair<size_t, size_t>, float*> mmat_b = read_from_file<float>(m_size, sizes[i], path_b, 0);
        float *etalon_lin_mat = read_etalon_from_file<float>(m_size, path_etalon);

        float *lin_seq_mult = seq_block_mat_multiplication<float>(mmat_a, mmat_b, m_size, sizes[i], time_seq);
        float *lin_int_mult = internal_parallel_block_mat_multiplication<float>(mmat_a, mmat_b, m_size, sizes[i], num_th, time_internal);
        float *lin_ext_mult = external_parallel_block_mat_multiplication<float>(mmat_a, mmat_b, m_size, sizes[i], num_th, time_external);

        bool compare_seq = compare_two_matrices<float>(lin_seq_mult, etalon_lin_mat, m_size, eps);
        bool compare_int = compare_two_matrices<float>(lin_int_mult, etalon_lin_mat, m_size, eps);
        bool compare_ext = compare_two_matrices<float>(lin_ext_mult, etalon_lin_mat, m_size, eps);

        if (!compare_seq)
            throw std::runtime_error("Wrong multiplication in 'seq_block_mat_multiplication function'\n");
        if (!compare_int)
            throw std::runtime_error("Wrong multiplication in 'internal_parallel_block_mat_multiplication function'\n");
        if (!compare_ext)
            throw std::runtime_error("Wrong multiplication in 'external_parallel_block_mat_multiplication function'\n");

        //push to lists
        res_seq.push_back(std::make_pair(sizes[i], time_seq));
        res_int.push_back(std::make_pair(sizes[i], time_internal));
        res_ext.push_back(std::make_pair(sizes[i], time_external));

        //print information
        std::cout << "Size block[ " << sizes[i] << " ][ float ]\n";
        std::cout << "Sequential multiplication time(s.):            " << time_seq << "\n";
        std::cout << "Internal parallel multiplication time(s.):     " << time_internal << "\n";
        std::cout << "External parallel multiplication time(s.):     " << time_external << "\n\n";

        release_map_matrix<float>(mmat_a, m_size, sizes[i], 1);
        release_map_matrix<float>(mmat_b, m_size, sizes[i], 0);
        release_linear_mtx<float>(etalon_lin_mat);

        release_linear_mtx<float>(lin_seq_mult);
        release_linear_mtx<float>(lin_int_mult);
        release_linear_mtx<float>(lin_ext_mult);
    }
    return std::make_tuple(res_seq, res_int, res_ext);
}

/**
 * @brief perform_double_test - test for double type
 * @param m_size
 * @param num_th
 * @return
 */
std::tuple<std::list<std::pair<size_t, double> >,
           std::list<std::pair<size_t, double> >,
           std::list<std::pair<size_t, double> > > perform_double_test(const size_t &m_size, const size_t num_th = 4)
{
    assert(m_size >= 2);

    std::list<std::pair<size_t, double> > res_seq, res_int, res_ext;
    for (size_t i = 1; i < sizes.size(); i++)
    {
        double time_seq = 0.0, time_internal = 0.0, time_external = 0.0;
        std::string path_a = "/home/rekrut/QTProjects/MultiplicationBlockMatrix/matrices/" + std::to_string(sizes[i]) + "/a.txt";
        std::string path_b = "/home/rekrut/QTProjects/MultiplicationBlockMatrix/matrices/" + std::to_string(sizes[i]) + "/b.txt";
        std::string path_etalon = "/home/rekrut/QTProjects/MultiplicationBlockMatrix/matrices/" + std::to_string(sizes[i]) + "/etalon.txt";
        std::map<std::pair<size_t, size_t>, double*> mmat_a = read_from_file<double>(m_size, sizes[i], path_a, 1);
        std::map<std::pair<size_t, size_t>, double*> mmat_b = read_from_file<double>(m_size, sizes[i], path_b, 0);
        double *etalon_lin_mat = read_etalon_from_file<double>(m_size, path_etalon);

        double *lin_seq_mult = seq_block_mat_multiplication<double>(mmat_a, mmat_b, m_size, sizes[i], time_seq);
        double *lin_int_mult = internal_parallel_block_mat_multiplication<double>(mmat_a, mmat_b, m_size, sizes[i], num_th, time_internal);
        double *lin_ext_mult = external_parallel_block_mat_multiplication<double>(mmat_a, mmat_b, m_size, sizes[i], num_th, time_external);

        bool compare_seq = compare_two_matrices<double>(lin_seq_mult, etalon_lin_mat, m_size, eps);
        bool compare_int = compare_two_matrices<double>(lin_int_mult, etalon_lin_mat, m_size, eps);
        bool compare_ext = compare_two_matrices<double>(lin_ext_mult, etalon_lin_mat, m_size, eps);

        if (!compare_seq)
            throw std::runtime_error("Wrong multiplication in 'seq_block_mat_multiplication function'\n");
        if (!compare_int)
            throw std::runtime_error("Wrong multiplication in 'internal_parallel_block_mat_multiplication function'\n");
        if (!compare_ext)
            throw std::runtime_error("Wrong multiplication in 'external_parallel_block_mat_multiplication function'\n");

        //push to lists
        res_seq.push_back(std::make_pair(sizes[i], time_seq));
        res_int.push_back(std::make_pair(sizes[i], time_internal));
        res_ext.push_back(std::make_pair(sizes[i], time_external));

        //print information
        std::cout << "Size block[ " << sizes[i] << " ][ double ]\n";
        std::cout << "Sequential multiplication time(s.):            " << time_seq << "\n";
        std::cout << "Internal parallel multiplication time(s.):     " << time_internal << "\n";
        std::cout << "External parallel multiplication time(s.):     " << time_external << "\n\n";

        release_map_matrix<double>(mmat_a, m_size, sizes[i], 1);
        release_map_matrix<double>(mmat_b, m_size, sizes[i], 0);

        release_linear_mtx<double>(etalon_lin_mat);
        release_linear_mtx<double>(lin_seq_mult);
        release_linear_mtx<double>(lin_int_mult);
        release_linear_mtx<double>(lin_ext_mult);
    }
    return std::make_tuple(res_seq, res_int, res_ext);
}

int main(int argc, char *argv[])
{
    std::srand( (unsigned) time(0) );

    //generate_files(m_size);

//    try
//    {
//        std::cout << "***************Start calculations\n";
//        std::tuple<std::list<std::pair<size_t, double> >,
//                   std::list<std::pair<size_t, double> >,
//                   std::list<std::pair<size_t, double> > > float_test = perform_float_test(m_size, 4);

//        std::tuple<std::list<std::pair<size_t, double> >,
//                   std::list<std::pair<size_t, double> >,
//                   std::list<std::pair<size_t, double> > > double_test = perform_double_test(m_size, 4);

//        create_csv_file(std::get<0>(float_test), "/home/rekrut/QTProjects/MultiplicationBlockMatrix/results/clang/float/seq_float.csv");
//        create_csv_file(std::get<1>(float_test), "/home/rekrut/QTProjects/MultiplicationBlockMatrix/results/clang/float/int_float.csv");
//        create_csv_file(std::get<2>(float_test), "/home/rekrut/QTProjects/MultiplicationBlockMatrix/results/clang/float/ext_float.csv");

//        create_csv_file(std::get<0>(double_test), "/home/rekrut/QTProjects/MultiplicationBlockMatrix/results/clang/double/seq_double.csv");
//        create_csv_file(std::get<1>(double_test), "/home/rekrut/QTProjects/MultiplicationBlockMatrix/results/clang/double/int_double.csv");
//        create_csv_file(std::get<2>(double_test), "/home/rekrut/QTProjects/MultiplicationBlockMatrix/results/clang/double/ext_double.csv");
//    }
//    catch (const std::runtime_error &re)
//    {
//        std::cerr << re.what() << "\n";
//    }
//    catch (...)
//    {
//        std::cout << "Unknown exception.\n";
//    }

    const size_t m_size = 2880;
    size_t b_size = 36;
    size_t num_th = 4;

    double time_seq = 0.0, time_internal = 0.0, time_external = 0.0;

    std::string path_a = "matrices/" + std::to_string(b_size) + "/a.txt";
    std::string path_b = "matrices/" + std::to_string(b_size) + "/b.txt";
    std::string path_etalon = "matrices/" + std::to_string(b_size) + "/etalon.txt";
    std::map<std::pair<size_t, size_t>, double*> mmat_a = read_from_file<double>(m_size, b_size, path_a, 1);
    std::map<std::pair<size_t, size_t>, double*> mmat_b = read_from_file<double>(m_size, b_size, path_b, 0);
    double *etalon_lin_mat = read_etalon_from_file<double>(m_size, path_etalon);

    double *lin_seq_mult = seq_block_mat_multiplication<double>(mmat_a, mmat_b, m_size, b_size, time_seq);
    double *lin_int_mult = internal_parallel_block_mat_multiplication<double>(mmat_a, mmat_b, m_size, b_size, num_th, time_internal);
    double *lin_ext_mult = external_parallel_block_mat_multiplication<double>(mmat_a, mmat_b, m_size, b_size, num_th, time_external);

    bool compare_seq = compare_two_matrices<double>(lin_seq_mult, etalon_lin_mat, m_size, eps);
    bool compare_int = compare_two_matrices<double>(lin_int_mult, etalon_lin_mat, m_size, eps);
    bool compare_ext = compare_two_matrices<double>(lin_ext_mult, etalon_lin_mat, m_size, eps);

    if (!compare_seq)
        throw std::runtime_error("Wrong multiplication in 'seq_block_mat_multiplication function'\n");
    if (!compare_int)
        throw std::runtime_error("Wrong multiplication in 'internal_parallel_block_mat_multiplication function'\n");
    if (!compare_ext)
        throw std::runtime_error("Wrong multiplication in 'external_parallel_block_mat_multiplication function'\n");

    //print information
    std::cout << "Size block[ " << b_size << " ][ double ]\n";
    std::cout << "Sequential multiplication time(s.):            " << time_seq << "\n";
    std::cout << "Internal parallel multiplication time(s.):     " << time_internal << "\n";
    std::cout << "External parallel multiplication time(s.):     " << time_external << "\n\n";

    return 0;
}
