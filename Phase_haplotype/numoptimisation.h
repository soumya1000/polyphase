#pragma once
#define NUMOPTIMISATION_H_INCLUDED

#include "globals.h"
#include "emalgo.h"



//methods for Numerical optimisation required by GSL 

  /* The gradient of f, df = (df/dx, df/dy). */
  double func_eval(const gsl_vector *v,  void *params);
  
  void diff_eval(const gsl_vector *v, void *params, gsl_vector *df);
    
  /* Compute both f and df together. */
  void func_diff_eval(const gsl_vector *x, void *params, double *f, gsl_vector *df);  
    
  int  optimise_param(em_hmm_genotype &);
  