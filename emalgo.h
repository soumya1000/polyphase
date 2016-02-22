#pragma once
#define EMALGO_H_INCLUDED
#include "globals.h"
#include "parsedata.h"
#include "particlefilter.h"
#include "sampling.h"
#include <cstdlib>
#include "haplo.h"
//to facilitate dakota communication
struct index_vars
{
  int marker;
  int cluster;
};

class em_hmm_genotype
{
    private:
   // ofstream m_params_file; 
    
    // This map is to track the gradient indices as interpreted by Dakota.
    vector<map<int,index_vars>> dakota_map;
    
    particle_filter m_part_filter;
    
    //list of the entire statespace mapping state-id to n-tuples of ancestral clusters
    //vector<states> m_total_states; //choose((N+K-1),N)
    vector<vector<int>> m_total_states; //choose((N+K-1),N)
    vector<vector<vector<int>>>  m_current_states; //  ind x markers x stateIDs
    vector<vector<vector<int>>>  m_current_transitions;
    //this is analogous to transition probs of a HMM  
    map<int,vector<Double>> m_jump_prob; // nMarkers x n_clusters (J=0,1 and 2)                  
  
    //this is analogous to emission probs of a HMM
    vector<vector<vector<Double>>> m_geno_given_clst_v;//nIndividuals x nMarkers x states
   
    //fwd probabilities
     vector<vector<vector<Double>>> m_fwd_probs;

    //Bckwd probabolities   
    vector<vector<vector<Double>>> m_bckwd_probs;
  
    //an intermediate variable
    vector<vector<vector<Double>>>  m_clst_given_geno_v;
    
    //number of jumps that occur to cluster k,J(imk)
    vector<vector<vector<Double>>> m_ev_jump_given_geno;//nIndividuals x nMarkers x no of Clusters
  
    vector<vector<Double>> m_clst_given_v; // nMarkers (x no states )  transistions
    //member functions
    
     //The first expression in Appendix C
    void compute_jump_prob(void);
  
    //Equation 10 & 11 in article  
    void compute_geno_given_clst_v();
  
    //Equation 8 and 9 in article
    void compute_clust_given_v();
    
    //Definition for p(g|v)in Appendix A of the article 
    Double compute_geno_given_v(int ind, int mark);
  
    //Appendix A
    void util_hmm_fwd(int marker,int state,int njumps,vector<vector<int>> &jumpData, 
			  vector<Double>&alpha_data);
    void util_hmm_bckwd(int marker,int state,int njumps,vector<vector<int>> &jumpData,
			vector<vector<Double>>&alpha_data);
    //void compute_fwd_prob_OLD();
    void compute_fwd_prob();
    //void compute_bckwd_prob_OLD();
      void compute_bckwd_prob();
  
    //Equation at the end of Appendix C
    void compute_clst_given_geno_v();
   //void util_get_clustOccurance(int cluster, int ntimes,vector<int>&states);
   //void compute_exp_jump_given_geno();
  //  void compute_jump_given_geno();
    
    void clear_variables();
    void update_HMM_param();
    void update_statespace();
   
    void log_param();
    
    Double  get_theta_gradient(int iCluster,int iCurMarker);
    Double  get_alpha_gradient(int iParamCount,int iCurMarker);
    Double  get_recomb_gradient(int iCurMarker);
    void updateConvergedState(); 
     
    void log_results(vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase);
    public:
    Chaplotypes m_hap_data;
    em_hmm_genotype();
    ~em_hmm_genotype();
  
    const vector<vector<vector<int>>>& getcurrentStates();
    void get_total_stateSpace();
  
    //generate randomly chosen states at each marker
    void get_current_stateSpace();
  
    
    //Numerical optimisation start    
    Double func_eval_local();    
   Double func_eval_local_OLD();    
   Double diff_eval_local(int index);  
   // void update_param(const gsl_vector *x,int status);
    //Numerical optimisation    End
    void initialise_model_param();
    void resolve_phase();
    void test_get_params();//testing purpose
};


    