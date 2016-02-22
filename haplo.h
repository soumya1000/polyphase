#pragma once
#define HAPLO_H_INCLUDED
#include "globals.h"
#include "parsedata.h"

class Chaplotypes
{
   
 private:
   void initialise_param(void); // initialises v(alpha,theta,r)
   
  
   
   public:
     Chaplotypes();
     ~Chaplotypes();
   void initial_call_actions(int num_clusters=8);
   void restore_prevstate();
   void print_param(void); 
   data_genotypes m_input;
   void log_param_hap();
   // nMarkers x nClusters x nClusters
   vector<map< pair<int,int>,Double>> m_trans_prob_bet_clst;
   
   //nClusters x nClusters
   vector<vector<Double>> m_alpha; 


   // nMarkers x nClusters 
   vector<vector<Double>> m_theta; 
   
   vector<Double> m_recombinations;
      
   
   //Equations 3 & 4 in the article
    void compute_trans_prob_bet_clst(void);
   
    void optimise_param_helper();
   
   //transitions between unordered lists of clusters
   Double compute_trans_prob_bet_clst_tuples(int frmMarker,vector<vector<int>>&fromstatePerms,vector<int>&toState);
 
};
/* To remove recombinations being estimated, comment recomb code in initial_call_acions and uncomment in final_call_actions. 
 take care of optimise_param_helper function and the pass_values function in main.cpp. in diff_eval_local comment call to get_recomb_gradient*/