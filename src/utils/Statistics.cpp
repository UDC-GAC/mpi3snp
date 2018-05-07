/*
 * This file is part of MPI3SNP.
 * Copyright (C) 2018 by Christian Ponte
 *
 * MPI3SNP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MPI3SNP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file Statistics.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Statistics class members implementation.
 */

#include "Statistics.h"
#include <mpi.h>

Statistics::Statistics() {
    pthread_mutex_init(&ints_mutex, nullptr);
    pthread_mutex_init(&longs_mutex, nullptr);
    pthread_mutex_init(&timers_mutex, nullptr);
}

Statistics::~Statistics() {
    pthread_mutex_destroy(&ints_mutex);
    pthread_mutex_destroy(&longs_mutex);
    pthread_mutex_destroy(&timers_mutex);
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

std::string Statistics::To_string() {
    std::string output("Statistics\n");
    for (auto item : timers){
        if (!std::get<2>(item)) {
            output += "\t" + std::get<0>(item) + ": " + std::to_string(std::get<1>(item)) + " seconds\n";
        }
    }

    for (auto item : ints){
        output += "\t" + item.first + ": " + std::to_string(item.second) + "\n";
    }
    for (auto item : longs){
        output += "\t" + item.first + ": " + std::to_string(item.second) + "\n";
    }

    return output;
}

std::vector<std::tuple<std::string, double, bool>>::const_iterator
Statistics::Find_timer_label(const std::string &label) {
    auto it = timers.begin();
    while (it < timers.end() && std::get<0>(*it).compare(label) != 0) {
        it++;
    }
    return it;
}