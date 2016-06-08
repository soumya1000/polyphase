#pragma once
#define HAPLO_H_INCLUDED
#include "globals.h"
#include "parsedata.h"

class Chaplotypes
{
  public:
    Chaplotypes();
   ~Chaplotypes();
   ofstream m_params_file_haplos;
   data_genotypes m_input;
   
   // nMarkers x nClusters x nClusters
   vector<map< pair<int,int>,Double>> m_trans_prob_bet_clst;
   
   //nClusters x nClusters
   vector<vector<double>> m_alpha; 
   
   // nMarkers x nClusters 
   vector<vector<double>> m_theta; 
   
   vector<double> m_recombinations;
   
   
   
   void initialise_param(string ipFileName); // initialises v(alpha,theta,r)
   
   void print_param(void);
   //Equations 3 & 4 in the article
   void compute_trans_prob_bet_clst(void);
   void log_param_hap();
   
   //transitions between unordered lists of clusters
   Double compute_trans_prob_bet_clst_tuples(int frmMarker,vector<vector<int>>&fromstatePerms,vector<int>&toState);
 
};
