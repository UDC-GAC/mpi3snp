//
// Created by christian on 11/10/17.
//

#ifndef MPI3SNP_ARG_PARSER_H
#define MPI3SNP_ARG_PARSER_H

#include <map>
#include <vector>
#include <thread>

class Arg_parser {
public:
    struct Arguments {
        std::string tped;
        std::string tfam;
        std::string output;
        std::vector<std::pair<unsigned int, unsigned int>> gpu_map;
        unsigned int cpu_threads = std::thread::hardware_concurrency();
        bool use_mi = true;
        unsigned int output_num = 10;
        bool debug;
        bool benchmarking;
    };

    class Finalize : public std::runtime_error {
    public:
        explicit Finalize(const std::string &message) : runtime_error(message) {};

        ~Finalize() override = default;
    };

    Arg_parser(int &argc, char **argv);

    Arguments get_arguments();

private:
    const int argc;
    char **argv;
};

#endif //MPI3SNP_ARG_PARSER_H
