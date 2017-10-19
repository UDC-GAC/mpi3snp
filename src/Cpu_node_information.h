//
// Created by christian on 18/10/17.
//

#ifndef MPI3SNP_CPU_NODE_INFORMATION_H
#define MPI3SNP_CPU_NODE_INFORMATION_H


#include "Node_information.h"

class Cpu_node_information : public Node_information {
public:
    Cpu_node_information();

    explicit Cpu_node_information(const void *ptr);

    std::string hostname() const override;

    std::string cpu() const override;

    long memory() const override;

    std::string mpi_library_version() const override;

    std::vector<int> processes() const override;

    std::vector<std::string> gpus() const override;

protected:
    size_t to_byteblock(void **ptr) const override;

    void add_processes(std::vector<int> processes) override;

private:
    static constexpr const char *cpu_name_bash_command = "cat /proc/cpuinfo | grep \"model name\" | uniq | "
            "cut -d':' -f2 | tail --bytes=+2";

    static constexpr const char *physical_memory_command = "cat /proc/meminfo | grep \"MemTotal\" | tr -d ' ' | "
            "cut -d':' -f2 | cut -d'k' -f1";

    static std::string get_cpu_info();

    static long get_physical_memory();

    std::string hostname_str;
    std::string cpu_str;
    long memory_size;
    std::string mpi_library;
    std::vector<int> process_list;
    std::vector<std::string> gpu_list;
};


#endif //MPI3SNP_CPU_NODE_INFORMATION_H
