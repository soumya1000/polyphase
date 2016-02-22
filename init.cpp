#include <cstdlib>
#include <iostream>
#include <fstream>
#include<regex>
#include <boost/multiprecision/mpfr.hpp>
#include "emalgo.h"
#include "globals.h"

//g++ -std=c++11 init.cpp emalgo.o particlefilter.o sampling.o parsedata.o haplo.o globals.o -o init -lgmpxx -lgmp  -lmpfr
void get_initial_vlaues(em_hmm_genotype &genoObj,Chaplotypes &haploObj);
void test_code();

int main(int argc, char *argv[])
{
  //mpfr_set_default_prec(prec);
  
  //double dThreshold = pow(10,-40);
  //double scaling_factor = pow(10,30);
  //test_code();
  Chaplotypes haploObj;
  int dScaleFactor=0;
  
  haploObj.initial_call_actions(); 
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

void test_code()
{
  Chaplotypes haploObj;
//  double dThreshold = pow(10,-40);
 // double scaling_factor = pow(10,30);
  int dScaleFactor=0;
  
  haploObj.initial_call_actions(); 
  std::remove("scaling.in");
  em_hmm_genotype HMM_geno_Obj;
  get_initial_vlaues(HMM_geno_Obj,haploObj);   

  Double obj_fun = HMM_geno_Obj.func_eval_local(); 
  
  
   /*std::ifstream myReadFile;
   string line;
   myReadFile.open("dakota_param.in");
   while(myReadFile.good())
   {
	getline(myReadFile,line); // get line from file
	   
	if(line.find("initial_point")!=string::npos) // search
	{
             break;
	}
    }

   line = line.substr( line.find_first_of(" "),line.length());
   line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
   int counter=0;

   string tempStr;
   vector<Double> params;
   while(counter < 83)
   {
	tempStr = line.substr( 0,line.find_first_of("\t"));
	params.emplace_back(stold(tempStr));
	line = line.substr( line.find_first_of("\t"),line.length());
	line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	++counter;
    }
   tempStr = line.substr( 0,line.find_first_of("\t"));
   params.emplace_back(stold(tempStr));

   haploObj.m_theta.clear();
   haploObj.m_alpha.clear();
   haploObj.m_recombinations.clear();
   counter=0;
     for(int iCount=0;iCount < n_markers; iCount++)
    {
      vector<Double> dTempVec;
       for(int jCount=0;jCount < n_clusters; jCount++)
       {
	  dTempVec.emplace_back(params[counter]);
	 // cout << params[counter] << " ";
	  ++counter;
      }
       haploObj.m_theta.emplace_back(dTempVec);
   }
     for(int iCount=0;iCount < n_markers; iCount++)
    {
      vector<Double> dTempVec;
       for(int jCount=0;jCount < n_clusters; jCount++)
       {
	  dTempVec.emplace_back(params[counter]);
	//  cout << params[counter] << " ";
	  ++counter;
      }
       haploObj.m_alpha.emplace_back(dTempVec);
   }
    for(int jCount=0;jCount < n_markers-1; jCount++)
    {
	  haploObj.m_recombinations.emplace_back(params[counter]);
	 // cout << params[counter] << " ";
	  ++counter;
     }
   obj_fun = HMM_geno_Obj.func_eval_local();  */
  Double gradients;  
  int num_deriv_vars = 2*n_markers*n_clusters + (n_markers-1);
  cout << endl <<"num_deriv_vars "<< num_deriv_vars << endl;
 for (int i=0; i<num_deriv_vars;++i)
 {
      gradients = HMM_geno_Obj.diff_eval_local(i);
      cout << gradients << ' ' ;
  }
  cout << endl; 
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

