#pragma once
#define GLOBALS_H_INCLUDED

#include <iostream>
#include <stdlib.h>  
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>  
#include <math.h> 
#include <fstream>
#include <random>
#include<iterator>
#include <numeric>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_multimin.h>
#include <boost/multiprecision/mpfr.hpp>
#include "tbb/tbb.h"
#include "tbb/parallel_for.h"
#include <tbb/concurrent_vector.h>
#include <tbb/combinable.h>
//#include <gmp.h>
#include<thread>
#define EPSILON    (1.0E-8)

using namespace std;
using namespace tbb;
namespace bmp = boost::multiprecision;
//definitions
//typedef long double Double;
typedef bmp::number<bmp::mpfr_float_backend<500>> Double;
//definitions

//#define num_states 10


//variable decalrations
extern int n_individuals;
extern int n_markers;
extern int n_clusters;
extern int n_ploidy;
extern int num_states;
extern int n_scaler;
extern Double dThreshold;
extern Double scaling_factor;
extern double Recomb_InitialValue;
extern int n_runs;
extern int n_iterations;
extern int n_transitions;
/*struct states
{
  vector<int> values;
  Double weight;
};*/

//function declarations

//void sample(vector<states> &input,  vector<int> &exclude ,vector<int> &output,int size_output);
//void sample(vector<states> &input, vector<int> &output,int size_output);

int factorial(int n);

long int binomialCoef(int N,int r);

Double binomialProb(int N,int r,Double p);

void vector_permutation(std::vector<int> &now,std::vector<int> next,
std::vector<vector<int>> &permutations);

//void gen_combi_states(vector< states > &T);
void gen_combi_states(vector< vector<int> > &T);

//void gen_weights_states(vector< states > &T );

long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector<vector<int>> &T);

int compare(const vector<int> &left, const vector<int> &right,vector<int> &out);

void vector_combination(vector<int > &input,vector<int> combination, int offset, int k,vector< vector<int>> &T);

bool small_vectors(const vector<int> &a,const vector<int> &b); 

//generates a map between two int vectors
void get_map_vectors(vector<vector<int>>&input1, vector<int> &input2, vector<vector<std::pair<int,int>>> &output);

Double  randZeroToOne();
//void print_values( int size, const gsl_vector *x);

bool inline pairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) 
{
      return firstElem.first < secondElem.first;
}

bool inline pairCompare1(const std::pair<int, Double>& firstElem, const std::pair<int, Double>& secondElem) 
{
      return firstElem.second > secondElem.second;
}

bool inline is_close(Double firstElem,  Double secondElem)
{
    return (secondElem/firstElem >0.1) ;
}

inline std::vector<int> erase_indices(const std::vector<int>& data, std::vector<size_t>& indicesToDelete)
{
    if(indicesToDelete.empty())
        return data;

    std::vector<int> ret;
    ret.reserve(data.size() - indicesToDelete.size());

    std::sort(indicesToDelete.begin(), indicesToDelete.end());

    // now we can assume there is at least 1 element to delete. copy blocks at a time.
    std::vector<int>::const_iterator iter_beginBlock = data.begin();
    for(std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++ it)
    {
        std::vector<int>::const_iterator iter_endBlock = data.begin() + *it;
        if(iter_beginBlock != iter_endBlock)
        {
            std::copy(iter_beginBlock, iter_endBlock, std::back_inserter(ret));
        }
        iter_beginBlock = iter_endBlock + 1;
    }

    // finally do not forget to copy last block.
    if(iter_beginBlock != data.end())
    {
        std::copy(iter_beginBlock, data.end(), std::back_inserter(ret));
    }

    return ret;
}