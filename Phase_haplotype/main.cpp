#include "globals.h"
#include "emalgo.h"
#include "numoptimisation.h"  

void test_function(em_hmm_genotype &HMM_geno_Obj);

int main(int argc, char *argv[])
{
	string ipFilename;
	if ( argc != 2 ) 
	{
	    // argv[0] is the program name
	    cout<<"usage: "<< argv[0] <<" <filename>\n";
	    exit(1);
	}
	else 
	{
	    // argv[1] is  input filename to open
	    ifstream the_file ( argv[1] );	    
	    if ( !the_file.is_open() )
	    {
	      cout<<"Could not open file\n";
	      exit(1);
	    }
	    else 
	    {
	      ipFilename =  argv[1];
	    }
	}
      int status;
      em_hmm_genotype HMM_geno_Obj(ipFilename);
      //test_function(HMM_geno_Obj);
      status = optimise_param(HMM_geno_Obj);
      HMM_geno_Obj.resolve_phase();      
      return 0;
}
 
void test_function(em_hmm_genotype &HMM_geno_Obj)
{
     Double objective_fuc;
     gsl_vector *x;
     int local_n_scaler=0;
     
     x = HMM_geno_Obj.optimise_param_helper();
     cout << "Params " << endl;
     print_values(x->size,x);
     gsl_vector *df = gsl_vector_alloc (x->size);
      
     objective_fuc = HMM_geno_Obj.func_eval_local(x);
     cout << endl<<"Gradients " << endl;
      HMM_geno_Obj.diff_eval_local(x,df);
     
     print_values(df->size,df);
      HMM_geno_Obj.resolve_phase();
}