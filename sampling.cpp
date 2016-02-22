#include "sampling.h"


Phase::Phase()
{
  
}
Phase::~Phase()
{
  
}

void Phase::resolvePhase(const vector<vector<vector<Double>>> &Fwd_probs,const vector<vector<vector<Double>>> &Clst_givenGenoV,
			  Chaplotypes &HmmObj,const vector<states> &totStateSpace,
		         const vector<vector<vector<int>>> &curStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)
{
    //  cout << "CHK Phase::resolvePhase START"<<endl;
      vector<int> chosenStates;
      vector<vector<Double>> Fwd_probs_Ind,Clst_givenGenoV_Ind;
      vector<vector<int>> stateSpace_Ind;
  
      for(int indCount=0;indCount< n_individuals;++indCount)
      {
	    vector<vector<int>> orderedStates_Ind;
	    vector<vector<int>> final_Phase_Ind;
	    
	    Fwd_probs_Ind = Fwd_probs[indCount];
	    Clst_givenGenoV_Ind = Clst_givenGenoV[indCount];
	    stateSpace_Ind = curStateSpace[indCount];
	    
	    sample_States(Fwd_probs_Ind,Clst_givenGenoV_Ind,stateSpace_Ind,HmmObj,totStateSpace,chosenStates);  
	    
	     resolve_StateOrder(HmmObj,totStateSpace,chosenStates,orderedStates_Ind);  
	    orderedStates.emplace_back(orderedStates_Ind);
	    
	    do_phasing(indCount,HmmObj,orderedStates_Ind,final_Phase_Ind);  
	    final_Phase.emplace_back(final_Phase_Ind);
	  //cout << endl;
	    chosenStates.clear();
      }
   //cout << "CHK Phase::resolvePhase END"<<endl;
}

void Phase::resolvePhase(const vector<vector<vector<Double>>> &Fwd_probs,const vector<vector<vector<Double>>> &Clst_givenGenoV,
			  Chaplotypes &HmmObj,const vector<states> &totStateSpace,
		         const vector<vector<int>> &curStateSpace, vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)
{
     // cout << "CHK Phase::resolvePhase START "<<endl;
      vector<int> chosenStates;
      vector<vector<Double>> Fwd_probs_Ind,Clst_givenGenoV_Ind;
      vector<vector<int>> stateSpace_Ind = curStateSpace; 
     
       for(int indCount=0;indCount< n_individuals;++indCount)
       {
	    //cout << indCount << endl;
	    vector<vector<int>> final_Phase_Ind;
	    vector<vector<int>> orderedStates_Ind;
	    Fwd_probs_Ind = Fwd_probs[indCount];
	    Clst_givenGenoV_Ind = Clst_givenGenoV[indCount];
	     sample_States(Fwd_probs_Ind,Clst_givenGenoV_Ind,stateSpace_Ind,HmmObj,totStateSpace,chosenStates); 
	     resolve_StateOrder(HmmObj,totStateSpace,chosenStates,orderedStates_Ind); 
	     orderedStates.emplace_back(orderedStates_Ind);
	     do_phasing(indCount,HmmObj,orderedStates_Ind,final_Phase_Ind);  
	     final_Phase.emplace_back(final_Phase_Ind);
	
	    chosenStates.clear();
      }
     // cout <<"CHK Phase::resolvePhase END " <<endl;
  
}
void Phase::sample_States(vector<vector<Double>> &Fwd_probs, vector<vector<Double>> &Clst_givenGenoV, 
		    vector<vector<int>> &curStateSpace, Chaplotypes &HmmObj,const vector<states> &totStateSpace,vector<int> &chosenStates)
{
     // cout << "CHK Phase::sample_States START"<<endl;
      vector<int> curMarkstates,fromStateVec,tempVec,toStateVec;
      vector<vector<int>> perms_frmState;
      int bestStateInd_loc,bestStateInd;
      vector<Double> curMarkValues;
      Double dfwdProb,dclstGV,dCurVal,dTemp;
      
      //first sample ziM ~ p(ziM|g,v)~p(g,ZiM,v) 
      bestStateInd_loc = distance(begin(Clst_givenGenoV[n_markers-1]), max_element(begin(Clst_givenGenoV[n_markers-1]), end(Clst_givenGenoV[n_markers-1]))); 
       bestStateInd = curStateSpace[n_markers-1][bestStateInd_loc];
      //cout << "bestStateInd " << bestStateInd << " ";
      chosenStates.emplace_back(bestStateInd);
      dCurVal = Clst_givenGenoV[n_markers-1][bestStateInd_loc];  
      
     // cout << n_markers-1 << " " << bestStateInd_loc << " " << bestStateInd <<endl;
      //choose zim recursively from markers M-1,...,1
      for(int markCount=n_markers-2; markCount>=0; --markCount)
      {
	      //we have value for previous marker in bestStateInd
	      fromStateVec = totStateSpace[bestStateInd].values;
	      curMarkstates = curStateSpace[markCount];
	      vector_permutation(fromStateVec,tempVec,perms_frmState);
		
	      //compute the value at this marker for all the states
	      for(int stateCount=0;stateCount< num_states;++stateCount)
	      {
		      dfwdProb = Fwd_probs[markCount][stateCount];
		      toStateVec = totStateSpace[curMarkstates[stateCount]].values;
		      
		      dclstGV= HmmObj.compute_trans_prob_bet_clst_tuples(markCount,perms_frmState,toStateVec);
		      
		      dTemp = dCurVal*dfwdProb;
		      dTemp *= dclstGV;
		      curMarkValues.emplace_back(dTemp);  
	      }
	    
	      //choose the state that has maximum value at this marker
	      bestStateInd_loc = distance(begin(curMarkValues), max_element(begin(curMarkValues), end(curMarkValues))); 
	      //we donot care for exact value of current Value, only the relative measure at a given marker
	      /*while(dCurVal < dThreshold)
	      {
		dCurVal *= scaling_factor;
	      }*/
	      dCurVal *= curMarkValues[bestStateInd_loc];      
	      bestStateInd = curMarkstates[bestStateInd_loc];
	     // cout << "bestStateInd " << bestStateInd << " ";
	      chosenStates.emplace(chosenStates.begin(),bestStateInd); 
	     // cout << markCount << " " << bestStateInd_loc << " " << bestStateInd <<endl;
	      curMarkValues.clear();
	      perms_frmState.clear();
	
      }
     
 /*  for(auto &kvp:chosenStates)
   {
     cout << kvp << " ";       
   }
   cout << endl;*/
 //cout << "CHK Phase::sample_States END"<<endl;
}

void Phase::resolve_StateOrder( Chaplotypes &HmmObj,const vector<states> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates)
{
      //get_best_orderToState(int frmMarker,vector<int>&frmState,vector<vector<int>>&tostatePerms,vector<int>&tobestState)
       cout << "CHK Phase::resolve_StateOrder START"<<endl;
      map<pair<int,int>,Double> markerData;  
      vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
      vector<vector<int>> tostatePerms;
      vector<Double> transValues;
      vector<std::pair<int,int>> tempState;
      std ::pair<int,int> correctPair;
      int position;
      Double tempTrans;
	
      vector<int> frmStateVec,toStateVec,tempVec;
	
      //retain the existing order for the first marker and so start from the second marker.   
      orderedStates.emplace_back(totStateSpace[chosenStates[0]].values);
      //cout << "markCount 0 " << chosenStates[0] <<endl;
      
      //from second marker to the last marker, recursively determine the order
      for(int markCount=1;markCount<n_markers;++markCount)
      {
	      vector<int> curMarkstates;
	      
	      frmStateVec = orderedStates[markCount-1];
	      toStateVec =  totStateSpace[chosenStates[markCount]].values; 
	    
	      markerData = HmmObj.m_trans_prob_bet_clst[markCount-1];
	      //cout << "markCount " << markCount << " " <<chosenStates[markCount]<<" from state" << endl;
	      /* for(auto &k:frmStateVec)
		  cout << k << ",";
	      cout << " ** ";
	      for(auto &k:toStateVec)
		  cout << k << ",";
		cout << endl;*/
	      vector_permutation(toStateVec,tempVec,tostatePerms);
	      get_map_vectors(tostatePerms, frmStateVec, Map_Valid_stateTuples);
	    
	      for(auto &kvp:Map_Valid_stateTuples)
	      {
		      tempTrans = 1.0;
		      for(auto kv:kvp)
		      {
			  correctPair.first = kv.second;
			  correctPair.second = kv.first;
			  //cout << kv.first << "," << kv.second << " ";
			  tempTrans = exp(log(tempTrans) + log(markerData[correctPair]));
		      }
		      transValues.emplace_back(tempTrans);
		      //cout << "  ***  " << tempTrans << endl;
	      }
	  
	      position =distance(begin(transValues), max_element(begin(transValues), end(transValues)));
	      tempState = Map_Valid_stateTuples[position];
	    /* cout << " position" <<position << " ";
	      cout << "state " ;*/
	      for(auto &kv:tempState)
	      {
		    curMarkstates.emplace_back(kv.first);
		//cout << kv.first << ",";
	      }
	      //cout << endl;
	      orderedStates.emplace_back(curMarkstates);
	      Map_Valid_stateTuples.clear();
	      transValues.clear();
	      tostatePerms.clear();
      }
      
    /*for(auto &kvp:orderedStates)
      {
	for(auto kv:kvp)
	  cout << kv << ",";     
	  cout << "  ";
      }
    */  
   // cout << "CHK Phase::resolve_StateOrder END"<<endl;
}
  
void Phase::do_phasing(int indCount, Chaplotypes &HmmObj, 
vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase)
{
     //  cout << "CHK Phase::do_phasing START"<<endl;
      vector<vector<int>> genotype_Ind,permuted_Chaplotypes;  
      vector<int> geno_MarkerVec,cur_state,tempVec;
      vector< vector<std::pair<int,int>>> Map_Haplos_states;
      vector<Double> probValues,thetaCurMarker;
      Double perm_values,dTemphaplo,dThetaOne;
      vector<std::pair<int,int>> tempState;
      
      //long  double scaling_factor = pow(10,20);
      //long  double dThreshold = pow(10,-20);
      int iscaledTimes,ihaplo,position;
	
      genotype_Ind= HmmObj.m_input.m_genotypes[indCount];
      
      for(int markCount=0;markCount< n_markers;++markCount)
      {
	      vector<int> phased_MarkerVec;
	      geno_MarkerVec = genotype_Ind[markCount];
	      thetaCurMarker = HmmObj.m_theta[markCount];
	      vector_permutation(geno_MarkerVec,tempVec,permuted_Chaplotypes);
	      cur_state = orderedStates[markCount];       
	      // map the state and haplo combi in unique way, for eg haplo={0, 0, 1, 1} state={2,3,3,4} map= {(0, 2),(0, 3),(1, 3),(1, 4)}
	      get_map_vectors(permuted_Chaplotypes, cur_state, Map_Haplos_states);
	      
	      for(auto &kvp:Map_Haplos_states)
	      {
		  perm_values=1.0;
		  iscaledTimes =0;
		  
		  for(auto clst_it:kvp)
		  {
		      dThetaOne = thetaCurMarker[clst_it.second];
		      ihaplo =   clst_it.first;      
		      // Equation (1) in Scheet & Stephans 2006  article
		      //dTemphaplo = pow(dThetaOne,ihaplo) * pow((1.0-dThetaOne),(1.0-ihaplo));  
		      dTemphaplo = pow(dThetaOne,ihaplo);
		      dThetaOne = 1.0-dThetaOne;
		      ihaplo=1.0-ihaplo;  
		      dTemphaplo *= pow(dThetaOne,ihaplo);
		      perm_values = perm_values*dTemphaplo; 
		  }
		
		    //collect values from all possible permutations  
		    probValues.emplace_back(perm_values);	      
	  }
	  position =distance(begin(probValues), max_element(begin(probValues), end(probValues)));
	  tempState = Map_Haplos_states[position];
      
	  for(auto &kv:tempState)
	  {
		phased_MarkerVec.emplace_back(kv.first);
	  }
	  
	  final_Phase.emplace_back(phased_MarkerVec);
	  permuted_Chaplotypes.clear();
	  Map_Haplos_states.clear();
	  probValues.clear();     
      }
    /* for(auto &kvp:final_Phase)
      {
	for(auto kv:kvp)
	  cout << kv << ",";     
	  cout << "  ";
      }
      cout << endl;
      */
  //  cout << "CHK Phase::do_phasing END"<<endl;
}