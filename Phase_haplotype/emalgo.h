#pragma once
#define EMALGO_H_INCLUDED
#include "globals.h"
#include "haplo.h"
#include "parsedata.h"
#include "particlefilter.h"
#include "sampling.h"
#include <cstdlib>


class em_hmm_genotype
{
    private:
    int m_func_Scaler; // variable to keep track of scaling in objective function
    ofstream m_params_file;  
    
    Chaplotypes m_hap_data;
    
    particle_filter m_part_filter;
    
    //list of the entire statespace mapping state-id to n-tuples of ancestral clusters
    vector<states> m_total_states; //choose((N+K-1),N)
    vector<Double> m_stateweights;
    vector<vector<vector<int>>>  m_current_states; //  ind x markers x stateIDs
  
    //this is analogous to transition probs of a HMM  
    map<int,vector<Double>> m_jump_prob; // nMarkers x n_clusters (J=0,1 and 2)                  
  
    //this is analogous to emission probs of a HMM
    vector<vector<vector<Double>>> m_geno_given_clst_v;//nIndividuals x nMarkers x states
   
    //fwd probabilities
    vector<vector<vector<Double>>> m_fwd_probs;

    //Bckwd probabolities   
    vector<vector<vector<Double>>> m_bckwd_probs;
  
    //an intermediate variable
    vector<vector<vector<Double>>>  m_clst_given_geno_v;//nIndividuals x nMarkers x no states
   
    //number of jumps that occur to cluster k,J(imk)
   // vector<vector<vector<double>>> m_ev_jump_given_geno;//nIndividuals x nMarkers x no of Clusters
  
    vector<vector<vector<Double>>> m_clst_given_v; //nIndividuals x nMarkers x no states
    
    //member functions
    
     //The first expression in Appendix C
    void compute_jump_prob(void);
  
    //Equation 10 & 11 in article  
    void compute_geno_given_clst_v(void);
  
    //Equation 8 and 9 in article
    void compute_clust_given_v();
    
    //Definition for p(g|v)in Appendix A of the article 
    Double compute_geno_given_v(int ind, int mark);
  
    //Appendix A
    void util_hmm_fwd(int individual, int marker,int state,int njumps,vector<vector<int>> &jumpData, 
			  vector<Double>&alpha_data);
    void util_hmm_bckwd(int individual,int marker,int state,int njumps,vector<vector<int>> &jumpData,
			vector<vector<Double>>&alpha_data);
    void compute_fwd_prob(void);
    void compute_bckwd_prob(void);
  
    //Equation at the end of Appendix C
    void compute_clst_given_geno_v(void);
    void compute_jump_given_geno(void);
    
    void clear_variables();
    void update_HMM_param();
       
    void get_theta_gradient(gsl_vector *df,int *iParamCount);//,int iCurMarker);
    void get_alpha_gradient(gsl_vector *df,int *iParamCount);//,int iCurMarker);
    void get_recomb_gradient(gsl_vector *df,int *iParamCount);//,int iCurMarker);
    
    void log_results(vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase);
    public:
    
    em_hmm_genotype(string ipFileName);
    ~em_hmm_genotype();
  
    gsl_vector* optimise_param_helper();
    void update_statespace(const gsl_vector *v);
    void log_param();
    void get_total_stateSpace();
  
    //generate randomly chosen states at each marker
    void get_current_stateSpace();  
  
    
    //Numerical optimisation start    
    double func_eval_local(const gsl_vector *v);    
    void diff_eval_local(const gsl_vector *v,gsl_vector *df);  
    void update_param(const gsl_vector *x,int status);
    //Numerical optimisation    End
    void initialise_model_param(string ipFileName);
    void resolve_phase();
      
};


    