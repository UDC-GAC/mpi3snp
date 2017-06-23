//
// Created by christian on 09/06/17.
//

#ifndef MPI3SNP_STATISTICS_H
#define MPI3SNP_STATISTICS_H

#include <vector>
#include <tuple>
#include <string>
#include <pthread.h>

class Statistics {
public:
    Statistics();

    ~Statistics();

    void Begin_timer(const std::string &label);

    double End_timer(const std::string &label);

    double Get_timer(const std::string &label);

    std::vector<std::pair<std::string, double>> Get_all_timers();

    void Add_value(const std::string &label, int value);

    int Get_value(const std::string &label);

    std::vector<std::pair<std::string, int>> Get_all_values();

    std::string To_string();

private:
    std::vector<std::tuple<std::string, double, bool>>::const_iterator Find_timer_label(const std::string &label);

    std::vector<std::pair<std::string, int>>::const_iterator Find_value_label(const std::string &label);

    pthread_mutex_t timers_mutex, values_mutex;
    std::vector<std::tuple<std::string, double, bool>> timers;
    std::vector<std::pair<std::string, int>> values;
};

#endif //MPI3SNP_STATISTICS_H
