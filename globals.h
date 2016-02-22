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
#include <boost/multiprecision/mpfr.hpp>
#include "tbb/tbb.h"
#include "tbb/parallel_for.h"
#include <tbb/concurrent_vector.h>
#include <tbb/combinable.h>
//#include <chrono>
//#include "./gmpfrxx/gmpfrxx.h"
#define EPSILON    (1.0E-8)
using namespace std;
using namespace tbb;
namespace bmp = boost::multiprecision;
//definitions
//typedef long double Double;
typedef bmp::number<bmp::mpfr_float_backend<500>> Double;
//typedef  bmp::mpfr_float  Double;
const Double Recomb_InitialValue = 0.06;
const Double Threshold_transitions= pow(10,-30);
const int prec =1000;


//variable decalrations
extern int n_individuals;
extern int n_markers;
extern int n_clusters;
extern int n_ploidy;
extern   Double dThreshold;
extern   Double scaling_factor;
extern int num_states;
struct states
{
  vector<int> values;
  Double weight;
};

//function declarations

void sample(vector<states> &input,  vector<int> &exclude ,vector<int> &output,int size_output);
void sample(vector<states> &input, vector<int> &output,int size_output);

int factorial(int n);

long int binomialCoef(int N,int r);

Double binomialProb(int N,int r,Double p);

void vector_permutation(std::vector<int> &now,
			
std::vector<int> next,std::vector<vector<int>> &permutations);

//void gen_combi_states(int nploidy,int nclusters,vector< states > &T);
void gen_combi_states(int nploidy,int nclusters,vector< vector<int> > &T);
void gen_weights_states(int n, int k, vector< states > &T ,vector<Double> &indiprobs);

long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector<vector<int>> &T);
long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector< vector<int>> &T);
int compare(const vector<int> &left, const vector<int> &right,vector<int> &out);

void vector_combination(vector<int > &input,vector<int> combination, int offset, int k,vector< vector<int>> &T);

bool small_vectors(const vector<int> &a,const vector<int> &b); 

 
void normalize_vector(vector<double> &inout);

void project_range0_1(vector<double> &inout);

//generates a map between two int vectors
void get_map_vectors(vector<vector<int>>&input1, vector<int> &input2, vector<vector<std::pair<int,int>>> &output);
Double  randZeroToOne();

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
