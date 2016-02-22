#pragma once
#define SAMPLING_H_INCLUDED

#include "globals.h"
#include "emalgo.h"

class Chaplotypes;

class Phase
{
  
  /*void sample_States( vector<vector<Double>> &Fwd_probs, vector<vector<Double>> &Clst_givenGenoV, 
		   vector<vector<int>> &stateSpace, Chaplotypes &HmmObj,const vector<states> &totStateSpace,vector<int> &chosenStates);*/
  void sample_States_OLD(vector<vector<Double>> &Fwd_probs, vector<vector<Double>> &Clst_givenGenoV, 
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates);
  void resolve_StateOrder_OLD( Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates);
   void do_phasing_OLD(int indCount, Chaplotypes &HmmObj, vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase);
  void sample_States( vector<vector<Double>> &Fwd_probs, vector<vector<Double>> &Clst_givenGenoV, 
		 Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates);
  void resolve_StateOrder( Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates);  
  void do_phasing(int indCount, Chaplotypes &HmmObj, vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase);


  public:
    
  Phase();
  ~Phase();

  /*  void resolvePhase(const vector<vector<vector<Double>>> &Fwd_probs,const vector<vector<vector<Double>>> &Clst_givenGenoV,
		    Chaplotypes &HmmObj,const vector<states> &totStateSpace,const vector<vector<vector<int>>> &curStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase);*/
  void resolvePhase(const vector<vector<vector<Double>>> &Fwd_probs,const vector<vector<vector<Double>>> &Clst_givenGenoV,
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase);
};
