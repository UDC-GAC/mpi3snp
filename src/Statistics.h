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
 * @file Statistics.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Statistics class definition and private template functions implementation.
 */

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

    inline void Addi(const std::string &label, int value) {
        pthread_mutex_lock(&ints_mutex);
        Add<int>(label, value, ints);
        pthread_mutex_unlock(&ints_mutex);
    }

    inline int Geti(const std::string &label) {
        int retval = 0;
        pthread_mutex_lock(&ints_mutex);
        retval = Get<int>(label, ints);
        pthread_mutex_unlock(&ints_mutex);
        return retval;
    }

    inline void Addl(const std::string &label, long value) {
        pthread_mutex_lock(&longs_mutex);
        Add<long>(label, value, longs);
        pthread_mutex_unlock(&longs_mutex);
    }

    inline long Getl(const std::string &label) {
        long retval = 0;
        pthread_mutex_lock(&longs_mutex);
        retval = Get<long>(label, longs);
        pthread_mutex_unlock(&longs_mutex);
        return retval;
    }

    void Begin_timer(const std::string &label);

    double End_timer(const std::string &label);

    double Get_timer(const std::string &label);

    std::string To_string();

private:
    template<typename T>
    static typename T::const_iterator Find(const std::string &label, const T &vector) {
        auto it = vector.begin();
        while (it < vector.end() && std::get<0>(*it).compare(label) != 0) {
            it++;
        }
        return it;
    };

    template<typename T>
    static void Add(const std::string &label, const T &value, std::vector<std::pair<std::string, T>> &vector) {
        auto pos = Find<std::vector<std::pair<std::string, T>>>(label, vector);
        if (pos != vector.end()) {
            pos = vector.erase(pos);
        }
        vector.insert(pos, std::make_pair(label, value));
    };

    template<typename T>
    static T Get(const std::string &label, const std::vector<std::pair<std::string, T>> &vector) {
        auto pos = Find<std::vector<std::pair<std::string, T>>>(label, vector);
        if (pos != vector.end()) {
            return std::get<1>(*pos);
        }
        return 0;
    };

    pthread_mutex_t ints_mutex;
    std::vector<std::pair<std::string, int>> ints;

    pthread_mutex_t longs_mutex;
    std::vector<std::pair<std::string, long>> longs;

    std::vector<std::tuple<std::string, double, bool>>::const_iterator Find_timer_label(const std::string &label);

    pthread_mutex_t timers_mutex;
    std::vector<std::tuple<std::string, double, bool>> timers;
};

#endif //MPI3SNP_STATISTICS_H
