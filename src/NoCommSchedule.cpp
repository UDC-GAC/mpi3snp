//
// Created by christian on 16/03/17.
//

#include "NoCommSchedule.h"

void NoCommSchedule::distribute_problem(size_t problem_size, unsigned int num_proc) {
    size_t block_size = problem_size / num_proc;
    int remainder = problem_size % num_proc;

    auto *pair_list = &dist[0];
    pair_list->reserve(num_proc);
    for (int i = 0; i < num_proc; i++) {
        std::pair<std::multiset<Block, std::less<Block>>, unsigned int> pair;
        size_t beginning = i > 0 ? (*pair_list)[i - 1].first.begin()->x + (*pair_list)[i - 1].first.begin()->xlen : 0;
        size_t size = block_size + (remainder-- > 0 ? 1 : 0);
        Block b(beginning, size);
        pair.first.insert(b);
        pair.second = i;
        pair_list->push_back(pair);
    }
}

void NoCommSchedule::combine_blocks(unsigned int k) {
    // Generate all combinations from the initial problem distribution
    for (int i = 1; i < k; i++) {
        dist[i].reserve(100);
        // For each (k-1)-combination
        for (auto c = dist[i - 1].begin(); c != dist[i - 1].end(); c++) {
            // Prefix it with the appropriate base element (which are all the base elements less than or equal to the
            // leftmost component)
            std::multiset<Block, std::less<Block>> set = c->first;
            Block first = *set.begin();
            for (auto block_iterator = dist[0].begin();
                 block_iterator != dist[0].end() && *block_iterator->first.begin() <= first;
                 block_iterator++) {
                std::pair<std::multiset<Block, std::less<Block>>, unsigned int> p;
                p.first.insert(*block_iterator->first.begin());
                p.first.insert(set.begin(), set.end());
                dist[i].push_back(p);
            }
        }
    }
}

void NoCommSchedule::assign_combinations() {

}

NoCommSchedule::NoCommSchedule(size_t problem_size, unsigned int k, unsigned int num_proc) {
    // Allocate distribution array
    dist = new std::vector<std::pair<std::multiset<Block, std::less<Block>>, unsigned int>>[k];

    // Divide the initial problem in blocks of contiguous data
    distribute_problem(problem_size, num_proc);
    // Generate the k-combinations of the different blocks
    combine_blocks(k);
    // Assign the combinations to the different processes
    assign_combinations();
}

NoCommSchedule::~NoCommSchedule() {}

int NoCommSchedule::Who_has(std::vector<Block *> comb) {
    return -1;
}

unsigned int NoCommSchedule::Who_has2(std::multiset<Block, std::less<Block>> c) {
    if (c.size() > k) {
        // Throw exception
        return 100;
    }

    auto selection = dist[c.size() - 1];

    // Binary search
    unsigned int lower = 0,
            upper = selection.size() - 1,
            index = (upper + lower) / 2;

    while (lower <= upper && selection[index].first != c) {
        if (Multiset_greater_than(c, selection[index].first)) {
            lower = index + 1;
        } else {
            upper = index - 1;
        }
        index = (upper + lower) / 2;
    }

    if (lower > upper) {
        // Not found
        return 100;
    } else {
        return selection[index].second;
    }
}

/* Compare multisets m1 and m2 in reverse lexicographic order */
bool NoCommSchedule::Multiset_greater_than(std::multiset<Block> m1, std::multiset<Block> m2) {
    auto it1 = m1.rbegin(),
            it2 = m2.rbegin();

    // Advance both iterators until a difference is found
    while (*it1 == *it2) {
        it1++;
        it2++;
    };

    return *it1 > *it2;
}

std::vector<std::vector<Block *>> NoCommSchedule::Get_blocks(unsigned int proc_id) {
    std::vector<std::vector<Block *>> empty;
    return empty;
}