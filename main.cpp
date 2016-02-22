#include "globals.h"
#include "emalgo.h"
//#include "numoptimisation.h"  
#include <cctype>
#include <cstdlib>
#include <boost/multiprecision/cpp_dec_float.hpp>
  //Test code
 // dakota -i phase_param.in -o phase_param.out > phase_param_std.out
void pass_values(em_hmm_genotype &genoObj, map<int,Double> &var); 
  
int main(int argc, char *argv[])
{
  //mpfr_set_default_prec(prec);
      try
     {
	  em_hmm_genotype HMM_geno_Obj;
	  ifstream fin(argv[1]);
  
	//cout << "argv[1] " << argv[1] <<endl;
	if (!fin)
	{
	    cerr << "\nError: failure opening " << argv[1] << endl;
	    exit(-1);
         }
      size_t i, j, num_vars, num_fns, num_deriv_vars;
      string vars_text, fns_text, dvv_text;
 
  
    // Get the parameter vector and ignore the labels
  
      fin >> num_vars >> vars_text;
     // cout << "num_vars " << num_vars << " vars_text " << vars_text<<endl;
 
      map<int,Double> vars;
      vector<int> labels(num_vars);
      Double var_i; string label_i; 
      int v_i;
      map<string, int>::iterator v_iter;
  
      for (i=0; i<num_vars; i++) 
      {
	  fin >> var_i >> label_i;
	  transform(label_i.begin(), label_i.end(), label_i.begin(),(int(*)(int))tolower);
	  //cout << " var_i " << var_i << " label_i "<< label_i <<" ";
	  v_i = i;//v_iter->second;
	  vars[v_i] = var_i;
	  labels[i] = v_i;
      }
      
      //cout << endl;
      // Get the ASV vector and ignore the labels
      fin >> num_fns >> fns_text;
      //cout << "num_fns " << num_fns <<endl << "ASV[i]  " ;
      vector<short> ASV(num_fns);
  
      for (i=0; i<num_fns; i++)
      {
	  fin >> ASV[i];
	  //cout << ASV[i] << " " ;
	  fin.ignore(256, '\n');
      }
      //cout << endl;
      // Get the DVV vector and ignore the labels
      fin >> num_deriv_vars >> dvv_text;
      vector<int> DVV(num_deriv_vars);
     //cout << "num_deriv_vars " << num_deriv_vars <<endl;
      unsigned int dvv_i;
      
      for (i=0; i<num_deriv_vars; i++) 
      {
	  fin >> dvv_i;
	  fin.ignore(256, '\n');
	  DVV[i] = labels[dvv_i-1];
	//cout << "DVV[i] " << DVV[i] << endl;
      }
 
      vector<Double> x ;
      Double obj_func,gradients ;
      for ( i=0; i<num_vars; i++) 
      {
	  x.push_back(vars[i]);
      }
      pass_values(HMM_geno_Obj, vars);
      ofstream fout(argv[2]);
  
      if (!fout)
      {
	  cerr << "\nError: failure creating " << argv[2] << endl;
	  exit(-1);
      }
  
      fout.precision(6); // 16 total digits
      //fout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);
      fout.setf(ios::scientific);
      fout.setf(ios::right);
      fout.width(5);
   
      if (ASV[0] & 1) // **** f:
      {
	  obj_func = HMM_geno_Obj.func_eval_local();     
	  //cout << obj_func << endl;
	  fout << "                     " << obj_func  << " f\n";   
      }

      if (ASV[0] & 2)  // **** df/dx:
     {
	  fout << "[ ";
	  //cout << "gradients " << endl ;
	  for (i=0; i<num_deriv_vars; i++)
	  {
	      gradients = HMM_geno_Obj.diff_eval_local(DVV[i]);	
	      fout << gradients << ' ' ;
	     // cout << gradients << ' ' ;
	  }
	  fout << "]\n";
	  //cout << endl;
  }

   fout.flush();
   fout.close();
  }
   catch(const std::exception& e)
   {
	std::cout << e.what() << '\n';
   }
  return 0;
 }
  
void pass_values(em_hmm_genotype &genoObj, map<int,Double> &var)
{
  int ipointer=0;
  genoObj.m_hap_data.m_theta.clear();
  genoObj.m_hap_data.m_alpha.clear();
  genoObj.m_hap_data.m_recombinations.clear();
   
  //cout << "Theta values " << endl;
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      vector<Double> temp;
      for(int jCount=0; jCount < n_clusters ;++jCount,++ipointer)
      {
	  //cout << var[ipointer] << " " ;
	  temp.emplace_back(var[ipointer]);
      }
      genoObj.m_hap_data.m_theta.emplace_back(temp);
      // cout << endl;
   }
    
   //cout << "Alpha values " << endl;
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      vector<Double> temp;
      for(int jCount=0; jCount < n_clusters; ++jCount,++ipointer)
      {
	  //cout << var[ipointer] << " " ;
	   temp.emplace_back(var[ipointer]);   
      }
      genoObj.m_hap_data.m_alpha.emplace_back(temp);
      //cout << endl;
    }
  
   //cout << "Recomb values " << endl;
  for(int iCount=0;iCount < n_markers-1; iCount++)
   {
     //cout << var[ipointer] << " " ;
     genoObj.m_hap_data.m_recombinations.emplace_back(var[ipointer]);
     ++ipointer;	
   }
    //cout << endl; 
}