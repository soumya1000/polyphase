#pragma once
#define PARTICLEFILTER_H_INCLUDED
#include "globals.h"
#include "haplo.h"

class particle_filter
{
   
    //returns sum over all individuals at each marker
    void summary_Markervals(const vector<vector<vector<Double>>> &fwdVals,
	const vector<vector<vector<Double>>> &bckwdVals, vector<vector<Double >> &T);
  
  //retains the top one third states  and samples for new states surrounding these states after discarding states with low values
    void filter_states(vector<vector<int>> &sortedSts, const vector<states>&tot_Space,vector<vector<int>> &upd_StateSpace,Chaplotypes &HmmObj);
    
  //ranks the states sorted according the summary of values for all individuals at each markerstate
    void rank_Statescontrib(vector<vector<Double>> &Input,
			      const vector<vector<int>> &curStateSpace,
			      vector<vector<int>> &T);
  //retains the top three states w.r.t summary of both forward and backward values each and samples for new states 
  //surrounding these states after discarding states with low values
    void filter_states(vector<vector<int>> &sortedStsFwd,vector<vector<int>> &sortedStsBwd,
			 const vector<states> tot_Space,vector<vector<int>> &upd_StateSpace);
   Double  util_add_states(int stateId,int markCnt, const vector<states> tot_Space,vector<int> &prev_markSpace,Chaplotypes &HmmObj);
  
  public:
      
    particle_filter();       
    ~particle_filter();
    
    //when individuals share same state space
    
     void filter_run(const vector<vector<vector<Double>>> &Forward_values, 
		       const vector<vector<vector<Double>>> &Backward_values,
		       const vector<vector<int>> &curStateSpace,const vector< states > &tot_StateSpace,
		       vector<vector<int>> &upd_StateSpace,Chaplotypes &HmmObj);
     
     //when individuals have different state spaces
    void filter_run(const vector<vector<vector<Double>>> &Forward_values, 
		       const vector<vector<vector<Double>>> &Backward_values,
		       const vector<vector<vector<int>>> &curStateSpace,const vector< states > &tot_StateSpace,
		       vector<vector<vector<int>>> &upd_StateSpace,Chaplotypes &HmmObj);
     
      
};