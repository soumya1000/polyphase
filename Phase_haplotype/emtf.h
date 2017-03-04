#pragma once
#define EMTF_H_INCLUDED
#include "globals.h"
#include "haplo.h"
#include "parsedata.h"
#include "samplingtf.h"

class CEmtfgeno
{
    ofstream m_params_file; 
    
    Chaplotypes m_hap_data;
    vector<vector<int>> m_total_states; //choose((N+K-1),N)
     //this is analogous to transition probs of a HMM  
    map<int,vector<Double>> m_jump_prob; // nMarkers x n_clusters (J=0,1 and 2)       
    
    vector<vector<vector<int>>>  m_current_transitions;// (markers-1) X statesFrom X allstates
  
  //this is analogous to emission probs of a HMM
    vector<vector<vector<Double>>> m_geno_given_clst_v;//nIndividuals x nMarkers x states
   
    //fwd probabilities
    vector<vector<vector<Double>>> m_fwd_probs;

    //Bckwd probabolities   
    vector<vector<vector<Double>>> m_bckwd_probs;
  
    //an intermediate variable
    vector<vector<vector<Double>>>  m_clst_given_geno_v;//nIndividuals x nMarkers x no states
      
  
    //number of jumps that occur to cluster k,J(imk)
   vector<vector<vector<Double>>> m_ev_jump_given_geno;//nIndividuals x nMarkers x no of Clusters
  
    vector<Double>m_initial_state_prob;
  vector<vector<vector<Double>>> m_clst_given_v; //nMarkers x no statesx no states
    
    //member functions    
     
     //The first expression in Appendix C
    void compute_jump_prob(void);  
    //Equation 10 & 11 in article  
    void compute_geno_given_clst_v(void);  
    void compute_fwd_prob(void);
    void compute_bckwd_prob(void);
    //Equation 8 and 9 in article
    void compute_InitialstateProbs();
    void compute_clust_given_v();    
    
    void update_HMM_param();
    void compute_exp_jump_given_geno_v(void);
  
    //Appendix C
    void compute_theta_updated(vector<vector<double> > &T);
    void compute_alpha_updated(vector<vector<double> > &T);
    void compute_recombination_updated(vector<double> &T);
    bool test_likelihood(double old_likelihood,double new_likelihood);
    
      //Equation at the end of Appendix C
    void compute_clst_given_geno_v(void);
      
    void log_results(vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase);
    bool inline is_dclose(double a, double b, double epsilon = 0.001)
   {
    return abs(a - b) < epsilon;
   }
    double compute_expectation();
    void filter_transitions(int curIteration);
    void populate_transitions();
    void  util_hmm_fwd(int marker,int state,int njumps,vector<vector<int>> &jumpData, 
			  vector<Double>&alpha_data);
public:
    void log_param();
    void initialise_param();
    //Read simulated values as initial params
    void initialise_Simmodel_param(string paramsFileName,bool recombDirect);
  
    void log_model_param(string fName,bool recombDirect);
    double run_em(int iterations);
    
    CEmtfgeno(string ipFileName,string paramsFileName=NULL);
    ~CEmtfgeno();
    
      //Numerical optimisation start  
    
   
    //Numerical optimisation    End
    
    void resolve_phase(bool looging=false);
  
};