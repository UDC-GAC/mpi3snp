/*
 * Macros.h
 *
 *  Created on: 10/10/2013
 *      Author: jorge
 */

#ifndef MACROS_H_
#define MACROS_H_

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <zlib.h>
#include <math.h>
#include <algorithm>
#include <pthread.h>
#include <iterator>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

using namespace std;

#define VERSION "1.0"

// To determine the block size
// The number of pairs that will have all the 3-SNP interactions calculated in each block
#define NUM_PAIRS_BLOCK 5

// For reading the SNP
#define SNP_MAX_NAME_SIZE 128

// For initial size of the arrays
#define DEFAULT_NUM_INDS 1024
#define DEFAULT_NUM_SNPS 65536

// Genotype encoding
#define GENO_AA	0
#define GENO_Aa	1
#define GENO_aa	2

// If performing the 2-SNP analysis too
#define DOUBLE_SNP_ANALYSIS
#define TRIPLE_SNP_ANALYSIS

// Output
//#define DEBUG // To print the logs of debugging
#define BENCHMARKING // To measure the computational time for a paper

#endif /* MACROS_H_ */
