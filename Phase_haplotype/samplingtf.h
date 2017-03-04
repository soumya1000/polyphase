#pragma once
#define SAMPLINGTF_H_INCLUDED

#include "globals.h"
#include "haplo.h"

class CPhasetf
{
  
  void sample_States( vector<vector<Double>> &Fwd_probs, vector<vector<vector<Double>>>&Clust_givenV,vector<vector<vector<int>>> &transitions,
		 Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates);
  void resolve_StateOrder( Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates);  
  void do_phasing(int indCount, Chaplotypes &HmmObj, vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase);


  public:
    
  CPhasetf();
  ~CPhasetf();

  void resolvePhase(const vector<vector<vector<Double>>> &Fwd_probs,vector<vector<vector<Double>>> &Clust_givenV,vector<vector<vector<int>>> &transitions,
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase,bool logging=false);
};