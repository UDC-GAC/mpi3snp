//
// Created by christian on 18/10/17.
//

#ifndef MPI3SNP_NODE_INFORMATION_H
#define MPI3SNP_NODE_INFORMATION_H

#include <vector>
#include <string>

class Node_information {
public:
    class Builder {
    public:
        static Node_information *get_information();
    };

    Node_information();

    virtual std::string hostname() const =0;

    virtual std::string cpu() const =0;

    virtual long memory() const =0;

    virtual std::string mpi_library_version() const =0;

    virtual std::vector<int> processes() const =0;

    virtual std::vector<std::string> gpus() const =0;

    static std::vector<Node_information *> gather(int process);

protected:
    virtual size_t to_byteblock(void **ptr) const =0;

    virtual void add_processes(std::vector<int> processes) =0;

    static Node_information *build_from_byteblock(const void *ptr);

    std::string hardware_id;

private:
    static constexpr const char *hid_bash_command = "echo \"$(cat /sys/class/pci_bus/*/device/*/vendor)$(cat "
            "/sys/class/pci_bus/*/device/*/device)\" | md5sum | cut -d' ' -f 1";
};


#endif //MPI3SNP_NODE_INFORMATION_H
