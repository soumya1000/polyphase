#include "numoptimisation.h"


int optimise_param(em_hmm_genotype &HMM_geno_Obj)
{
  //set up the gsl variabes
   //cout<< "numoptimisation.cpp : optimise_param() START" << endl; 
   int iter = 0;
   int status, num_Params;
   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;  
   gsl_vector *x;
  
   gsl_multimin_function_fdf optimiser_func;
  
     
   x = HMM_geno_Obj.optimise_param_helper();
   num_Params = x->size;  
  
   optimiser_func.n =  num_Params;
   optimiser_func.f = func_eval;
   optimiser_func.df = diff_eval;
   optimiser_func.fdf = func_diff_eval;
   optimiser_func.params = &HMM_geno_Obj;
     
 
   T = gsl_multimin_fdfminimizer_conjugate_fr;//;
 
   s = gsl_multimin_fdfminimizer_alloc (T, num_Params);
  
   gsl_multimin_fdfminimizer_set (s, &optimiser_func, x, 0.01, 1e-100);  
  
   bool flag= false; 
 
   ofstream gsl_log_file("gsl_ouputlog.txt");  
   
  do{ 
       
	    ++iter;
	    status = gsl_multimin_fdfminimizer_iterate(s);     
	     if (status)
		  break;
	 
	    status = gsl_multimin_test_gradient (s->gradient, 1e-100);
	    
	    if (status == GSL_SUCCESS)
	    {
		gsl_log_file << "Maximum found :" << endl;  	    
	    }

	  gsl_log_file << iter<< " "<<gsl_strerror(status) << " "<< s->f << endl;
	  HMM_geno_Obj.update_statespace(s->x);
    }while(status == GSL_CONTINUE && iter <100);
   
    
  HMM_geno_Obj.update_param(s->x,0); //update with final params
 
  gsl_multimin_fdfminimizer_free(s);  
  gsl_vector_free (x);

 //cout<< "numoptimisation.cpp : optimise_param() END" << endl;
  return 0;
//cout<< "numoptimisation.cpp : optimise_param() END" << endl;
}
   
double func_eval(const gsl_vector *v,  void *params)
{
  //cout<< "numoptimisation.cpp : func_eval START" << endl;
  
  em_hmm_genotype *em_hmm_genotype_Obj = (em_hmm_genotype*)params;
  //cout<< "numoptimisation.cpp : func_eval END" << endl;
  return  em_hmm_genotype_Obj->func_eval_local(v);
}

/* The gradient of f, df = (df/dtheta_km, df/dalpha_km,df/dr_m). */
void diff_eval(const gsl_vector *v, void *params, gsl_vector *df)
{
  //cout<< "numoptimisation.cpp : diff_eval START" << endl;
  em_hmm_genotype *em_hmm_genotype_Obj = (em_hmm_genotype*)params;
  em_hmm_genotype_Obj->diff_eval_local(v,df);
  
 // cout<< "numoptimisation.cpp : diff_eval END" << endl;
}

void func_diff_eval(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
 //cout<< "numoptimisation.cpp : func_diff_eval START" << endl;
     *f = func_eval(x, params); 
    //cout << *f << endl;
    diff_eval(x, params, df);
 //cout<< "numoptimisation.cpp : func_diff_eval END" << endl;
}
