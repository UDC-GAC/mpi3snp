//
// Created by christian on 16/03/17.
//

#include "NoCommSchedule.h"

void NoCommSchedule::Create_blocks(size_t problem_size, unsigned int num_proc) {
    size_t block_size = problem_size / num_proc;
    int remainder = problem_size % num_proc;

    auto *pair_list = &dist[0];
    pair_list->reserve(num_proc);
    for (unsigned int i = 0; i < num_proc; i++) {
        std::pair<std::multiset<Block, std::less<Block>>, unsigned int> pair;
        size_t beginning = i > 0 ? (*pair_list)[i - 1].first.begin()->x + (*pair_list)[i - 1].first.begin()->xlen : 0;
        size_t size = block_size + (remainder-- > 0 ? 1 : 0);
        Block b(beginning, size);
        pair.first.insert(b);
        pair.second = INVALID_PROC;
        pair_list->push_back(pair);
    }
}

/* Factorial function */
unsigned long fact(unsigned long f){
    unsigned long prod = 1;
    for (unsigned long i = f; i > 0; i--){
        prod *= i;
    }
    return prod;
};

void NoCommSchedule::Combine_blocks(unsigned int k) {
    /* Generate all combinations from the initial problem distribution */
    for (unsigned int i = 1; i < k; i++) {
        /* Allocate the necessary space for all the possible k-combinations: C(b,k) = (k+b-1)!/k!(b-1)! */
        dist[i].reserve(fact(dist[0].size()+i)/(fact(i+1) * fact(dist[0].size()-1)));
        /* For each (k-1)-combination */
        for (auto c = dist[i - 1].begin(); c != dist[i - 1].end(); c++) {
            /* Suffix all combinations with the base elements equal or greater than the rightmost component */
            std::multiset<Block, std::less<Block>> set = c->first;
            Block last = *set.rbegin();
            auto block_iterator = dist[0].begin();
            while (*block_iterator->first.begin() != last) {
                block_iterator++;
            }
            for (; block_iterator != dist[0].end(); block_iterator++) {
                std::pair<std::multiset<Block, std::less<Block>>, unsigned int> p;
                p.first.insert(set.begin(), set.end());
                p.first.insert(p.first.end(), *block_iterator->first.begin());
                p.second = INVALID_PROC;
                dist[i].push_back(p);
            }
        }
    }
}

/* Translation table: for each combination of range k, create a pair of:
 * -> Id of the process who has the contingence table of the k-1 elements on the left of the multiset.
 * -> Id of the owner of the element on the rightmost of the multiset. */
void NoCommSchedule::Generate_translation_table(int range, std::vector<std::pair<unsigned int, unsigned int>> *output) {
    for (auto it = dist[range].begin(); it != dist[range].end(); it++) {
        /* Separate the last component from the rest */
        std::multiset<Block, std::less<Block>> remainder, rightmost;
        auto current = it->first.begin();
        auto previous = current++;
        while (current != it->first.end()) {
            remainder.insert(remainder.end(), *previous);
            previous = current++;
        }
        rightmost.insert(rightmost.end(), *previous);

        /* Find which process has each part */
        std::pair<unsigned int, unsigned int> p;
        p.first = Who_has(remainder);
        p.second = Who_has(rightmost);
        output->push_back(p);
    }
}

unsigned int NoCommSchedule::Optimal_selection(std::vector<std::pair<unsigned int, unsigned int>> *tr,
                                               unsigned int range, unsigned int p, int flag) {
    unsigned int i = 0;
    switch (flag) {
        /* Perfect match: the process has both the block and the contingence table */
        case 0:
            /* Advance until a unassigned perfect match is found */
            while (i < tr->size() &&
                   (dist[range][i].second != INVALID_PROC || !(tr->at(i).first == p && tr->at(i).second == p))) {
                i++;
            };
            break;

            /* Favourable match: the process has the contingence table */
        case 1:
            while (i < tr->size() && (dist[range][i].second != INVALID_PROC || tr->at(i).first != p)) {
                i++;
            };
            break;

            /* Unfavourable match: the process has the basic block */
        case 2:
            while (i < tr->size() && (dist[range][i].second != INVALID_PROC || tr->at(i).second != p)) {
                i++;
            };
            break;

            /* Mismatch: the process doesn't have neither the contingence table nor the basic block */
        default:
            while (i < tr->size() && dist[range][i].second != INVALID_PROC) {
                i++;
            };
            break;
    }

    if (i < tr->size()) {
        /* No combination matches the search */
        return i;
    } else {
        return INVALID_PROC;
    }
}

void NoCommSchedule::Assign_combinations(unsigned int num_proc) {
    /* Basic block distribution */
    unsigned int p;
    std::vector<std::pair<std::multiset<Block, std::less<Block>>, unsigned int>>::iterator it;
    for (it = dist[0].begin(), p = 0; it != dist[0].end(); it++, p = (p + 1) % num_proc) {
        it->second = p;
    }

    /* Block combination distribution */
    for (unsigned int i = 1; i < k; i++) {
        /* Generate translation table */
        std::vector<std::pair<unsigned int, unsigned int>> tr_table;
        Generate_translation_table(i, &tr_table);

        /* Distribute all combinations using the translation table */
        int proc_flag[num_proc];
        unsigned long count = dist[i].size();
        for (unsigned int j = 0; j < num_proc; j++)
            proc_flag[j] = 0;
        p = 0;
        while (count > 0) {
            unsigned int result;
            while ((result = Optimal_selection(&tr_table, i, p, proc_flag[p])) == INVALID_PROC) {
                proc_flag[p]++;
            }
            dist[i][result].second = p;
            p = (p + 1) % num_proc;
            count--;
        }
    }
}

NoCommSchedule::NoCommSchedule(size_t problem_size, unsigned int k, unsigned int num_proc) :
        k(k) {
    /* Allocate distribution array */
    dist = new std::vector<std::pair<std::multiset<Block, std::less<Block>>, unsigned int>>[k];

    /* Divide the initial problem in blocks of contiguous data */
    Create_blocks(problem_size, num_proc);
    /* Generate the k-combinations of the different blocks */
    Combine_blocks(k);
    /* Assign the combinations to the different processes */
    Assign_combinations(num_proc);
}

NoCommSchedule::~NoCommSchedule() {
    delete[] dist;
}

unsigned int NoCommSchedule::Who_has(std::multiset<Block, std::less<Block>> c) {
    if (c.size() > k) {
        /* TODO: return error via exception or code */
        return INVALID_PROC;
    }

    auto selection = dist[c.size() - 1];

    /* Binary search, since the combinations inside each vector are ordered lexicographically */
    unsigned long lower = 0,
            upper = selection.size() - 1,
            index = (upper + lower) / 2;

    while (lower <= upper && selection[index].first != c) {
        if (c > selection[index].first) {
            lower = index + 1;
        } else {
            upper = index - 1;
        }
        index = (upper + lower) / 2;
    }

    if (lower > upper) {
        return INVALID_PROC;
    } else {
        return selection[index].second;
    }
}

void NoCommSchedule::Get_blocks(unsigned int proc_id, std::vector<std::multiset<Block, std::less<Block>>> *block_list) {
    for (unsigned int i = 1; i < k; i++){
        for (auto it = dist[i].begin(); it != dist[i].end(); it++){
            if (it->second == proc_id){
                block_list->push_back(it->first);
            }
        }
    }
}