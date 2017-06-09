//
// Created by christian on 09/06/17.
//

#include "Statistics.h"
#include <mpi.h>

Statistics::Statistics() {
    pthread_mutex_init(&timers_mutex, NULL);
    pthread_mutex_init(&values_mutex, NULL);
}

Statistics::~Statistics() {
    pthread_mutex_destroy(&timers_mutex);
    pthread_mutex_destroy(&values_mutex);
}

std::vector<std::tuple<std::string, double, bool>>::const_iterator
Statistics::Find_timer_label(const std::string &label) {
    auto it = timers.begin();
    while (it < timers.end() && std::get<0>(*it).compare(label) != 0) {
        it++;
    }
    return it;
}

std::vector<std::pair<std::string, int>>::const_iterator Statistics::Find_value_label(const std::string &label) {
    auto it = values.begin();
    while (it < values.end() && it->first.compare(label) != 0) {
        it++;
    }
    return it;
}

void Statistics::Begin_timer(const std::string &label) {
    pthread_mutex_lock(&timers_mutex);
    auto pos = Find_timer_label(label);
    if (pos != timers.end()){
        pos = timers.erase(pos);
    }
    timers.insert(pos, std::make_tuple(label, MPI_Wtime(), true));
    pthread_mutex_unlock(&timers_mutex);
}

// TODO: raise exception if label does not exist
double Statistics::End_timer(const std::string &label) {
    pthread_mutex_lock(&timers_mutex);
    auto pos = Find_timer_label(label);
    if (pos != timers.end()) {
        double time;
        bool active;
        std::tie(std::ignore, time, active) = *pos;

        if (active) {
            time = MPI_Wtime() - time;
            pos = timers.erase(pos);
            timers.insert(pos, std::make_tuple(label, time, false));
            pthread_mutex_unlock(&timers_mutex);
            return time;
        }
    }
    pthread_mutex_unlock(&timers_mutex);
    return 0;
}

// TODO: raise exception if label does not exist
double Statistics::Get_timer(const std::string &label) {
    double out = 0;
    pthread_mutex_lock(&timers_mutex);
    auto pos = Find_timer_label(label);
    if (pos != timers.end()) {
        out = std::get<1>(*pos);
    }
    pthread_mutex_unlock(&timers_mutex);
    return out;
}

std::vector<std::pair<std::string, double>> Statistics::Get_all_timers() {
    std::vector<std::pair<std::string, double>> out;
    pthread_mutex_lock(&timers_mutex);
    for (auto it = timers.begin(); it < timers.end(); it++) {
        if (!std::get<2>(*it)) {
            out.push_back(std::make_pair(std::get<0>(*it), std::get<1>(*it)));
        }
    }
    pthread_mutex_unlock(&timers_mutex);
    return out;
}

void Statistics::Add_value(const std::string &label, int value) {
    pthread_mutex_lock(&values_mutex);
    auto pos = Find_value_label(label);
    if (pos != values.end()){
        pos = values.erase(pos);
    }
    values.insert(pos, std::make_pair(label, value));
    pthread_mutex_unlock(&values_mutex);
}

// TODO: raise exception if label does not exist
int Statistics::Get_value(const std::string &label) {
    int out = 0;
    pthread_mutex_lock(&values_mutex);
    auto pos = Find_value_label(label);
    if (pos != values.end()) {
        out = pos->second;
    }
    pthread_mutex_unlock(&values_mutex);
    return out;
}

std::vector<std::pair<std::string, int>> Statistics::Get_all_values() {
    std::vector<std::pair<std::string, int>> out;
    pthread_mutex_lock(&values_mutex);
    out = values;
    pthread_mutex_unlock(&values_mutex);
    return out;
}