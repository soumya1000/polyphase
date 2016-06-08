#pragma once
#define PARTICLEFILTER_H_INCLUDED
#include "globals.h"

class particle_filter
{
   
  //returns sum over all individuals at each marker
    //void summary_Markervals(const vector<vector<vector<Double>>> &Input,vector<vector<int>> &cur_states,vector<Double > &T);
      
  //ranks the states sorted according the summary of values for all individuals at each markerstate
    void rank_Statescontrib(vector<vector<Double>> &Input,
			      const vector<vector<int>> &curStateSpace,
			      vector<vector<int>> &T);
  //retains the top three states w.r.t summary of both forward and backward values each and samples for new states 
  //surrounding these states after discarding states with low values
    void filter_states(vector<vector<int>> &sortedSts, const vector<Double> tot_Space_weights,vector<vector<int>> &upd_StateSpace);
  
  public:
      
    particle_filter();       
    ~particle_filter();
    void filter_run(const vector<vector<vector<Double>>> &hmm_values, const vector<vector<vector<int>>> &curStateSpace,
		    const vector< Double > &tot_StateSpace_weights,vector<vector<vector<int>>> &upd_StateSpace);
    
      
};