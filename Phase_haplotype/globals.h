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
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/combinable.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
//#include <gmp.h>

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
struct states
{
  vector<int> values;
  Double weight;
};

//function declarations

void sample(vector<states> &input,  vector<int> &exclude ,vector<int> &output,int size_output);
void sample(vector<states> &input, vector<int> &output,int size_output);
void sample(vector<Double> &input, vector<int> &output,int size_output);
int factorial(int n);

long int binomialCoef(int N,int r);

Double binomialProb(int N,int r,Double p);

void vector_permutation(std::vector<int> &now,
			
std::vector<int> next,std::vector<vector<int>> &permutations);

void gen_combi_states(vector< states > &T);

void gen_weights_states(vector< states > &T );

long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector<vector<int>> &T);

int compare(const vector<int> &left, const vector<int> &right,vector<int> &out);

void vector_combination(vector<int > &input,vector<int> combination, int offset, int k,vector< vector<int>> &T);

bool small_vectors(const vector<int> &a,const vector<int> &b); 

//long double add_log(long double a,long double b);

//long double sub_log(long double a,long double b);

void normalize_vector(vector< Double> &inout);

//void project_range0_1(vector< Double> &inout);

//generates a map between two int vectors
void get_map_vectors(vector<vector<int>>&input1, vector<int> &input2, vector<vector<std::pair<int,int>>> &output);
Double  randZeroToOne();
void print_values( int size,  const gsl_vector *x);
