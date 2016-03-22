#pragma once
#define PARTICLEFILTER_H_INCLUDED
#include "globals.h"
#include "haplo.h"

class particle_filter
{
  
  private:
  //#########WHEN individuals share same state space and states are same for all markers ####
  //returns aggregate over all individuals and all markers
   void summary_Markervals(const vector<vector<vector<Double>>> &hmm_values, vector<Double> &T);
   
     //ranks the states sorted according the summary of values for all individuals 
    void rank_Statescontrib(vector<Double> &Input, const vector<int> &curStateSpace, vector<int> &T);
     //retains the top one third states  and samples for new states surrounding these states after discarding states with low values
    void filter_states(vector<int> &sortedSts, const vector<Double>&tot_Space_weights,vector<int> &upd_StateSpace,Chaplotypes &HmmObj);
    
    
  public:
      
    particle_filter();       
    ~particle_filter();
      //WHEN individuals share same state space and states are same for all markers
     void filter_run(const vector<vector<vector<Double>>> &hmm_values,
		       const vector<int> &curStateSpace,const vector< Double > &tot_StateSpace_weights,
		       vector<int> &upd_StateSpace,Chaplotypes &HmmObj);
     
         
};