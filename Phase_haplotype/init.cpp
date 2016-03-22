#include <cstdlib>
#include <iostream>
#include <fstream>
#include<regex>
#include <boost/multiprecision/mpfr.hpp>
#include "emalgo.h"
#include "globals.h"

//g++ -std=c++11 init.cpp -I/run/media/root/System/tbb44_20150728oss/include/ emalgo.o particlefilter.o sampling.o parsedata.o haplo.o globals.o -o init -lmpfr -L /run/media/root/System/tbb44_20150728oss/build/linux_intel64_gcc_cc4.8_libc2.19_kernel3.16.7_release/ -ltbb

void get_initial_vlaues(em_hmm_genotype &genoObj,Chaplotypes &haploObj);

int main(int argc, char *argv[])
{
 
  string ipFilename;
  
   if ( argc != 2 )  
   {
	cout<<"usage: "<< argv[0] <<" <filename>\n";
	exit(1);
   }
    else 
    {
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
 
  Chaplotypes haploObj;
  int dScaleFactor=0;
  
  haploObj.initial_call_actions(ipFilename); 
  std::remove("scaling.in");
  std::remove("poly_phase.dat");
  em_hmm_genotype HMM_geno_Obj;
  get_initial_vlaues(HMM_geno_Obj,haploObj);   
  
  Double obj_fun = HMM_geno_Obj.func_eval_local(); 
  Double obj_fun_copy = obj_fun;
  
  while(abs(obj_fun)<dThreshold)
  {
     obj_fun *= scaling_factor;
     ++dScaleFactor ;
  }
 
  std::ofstream scale_file("scaling.in");
  scale_file << dScaleFactor << endl;
  scale_file << obj_fun_copy << endl;
  scale_file << obj_fun << endl;
  scale_file << 1 << endl;
  scale_file.close();
  
  if(system("dakota -i dakota_param.in -o dakota_param.out > dakota_param_std.out") == 0) 
  {
	HMM_geno_Obj.update_statespace(); 
        HMM_geno_Obj.resolve_phase();
  }
  else
  {
      cout << "FAILURE: PARAMETER OPTIMIZATION" << endl;
  }
  return 0;
}

void get_initial_vlaues(em_hmm_genotype &genoObj,Chaplotypes &haploObj)
{
  
    genoObj.m_hap_data.m_theta.clear();
    genoObj.m_hap_data.m_alpha.clear();
    genoObj.m_hap_data.m_recombinations.clear();
  
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      genoObj.m_hap_data.m_theta.emplace_back(haploObj.m_theta[iCount]);
   }
    
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      genoObj.m_hap_data.m_alpha.emplace_back(haploObj.m_alpha[iCount]);
   }

  for(int iCount=0;iCount < n_markers-1; iCount++)
  {
      genoObj.m_hap_data.m_recombinations.emplace_back(haploObj.m_recombinations[iCount]);
  }

}

