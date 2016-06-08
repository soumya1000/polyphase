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
  cout << "CHK: Phase::resolvePhase() START" << endl;
  vector<int> chosenStates;
  vector<vector<Double>> Fwd_probs_Ind,Clst_givenGenoV_Ind;
  vector<vector<int>> stateSpace_Ind;
  
  try
  {
	for(int indCount=0;indCount< n_individuals;++indCount)
	{
	    vector<vector<int>> orderedStates_Ind;
	    vector<vector<int>> final_Phase_Ind;
	    
	    Fwd_probs_Ind = Fwd_probs[indCount];
	    Clst_givenGenoV_Ind = Clst_givenGenoV[indCount];
	    stateSpace_Ind = curStateSpace[indCount];
	    
	    //cout << "sample_States "<< indCount << endl;
	    
	    sample_States(Fwd_probs_Ind,Clst_givenGenoV_Ind,stateSpace_Ind,HmmObj,totStateSpace,chosenStates);  
	    
	    //cout << "resolve_StateOrder "<< indCount << endl;
	    resolve_StateOrder(HmmObj,totStateSpace,chosenStates,orderedStates_Ind);  
	    orderedStates.emplace_back(orderedStates_Ind);
	    
	    //cout << "do_phasing "<< indCount << endl;
	    do_phasing(indCount,HmmObj,orderedStates_Ind,final_Phase_Ind);  
	    final_Phase.emplace_back(final_Phase_Ind);
	  //cout << endl;
	    chosenStates.clear();
	}
  }  
  
  catch(const std::exception& e)
  {
	std::cout << e.what() << '\n';
  }
  cout << "CHK: Phase::resolvePhase() END" << endl;
}

void Phase::sample_States(vector<vector<Double>> &Fwd_probs, vector<vector<Double>> &Clst_givenGenoV, 
		    vector<vector<int>> &curStateSpace, Chaplotypes &HmmObj,const vector<states> &totStateSpace,vector<int> &chosenStates)
{
	//cout << "CHK: Phase::sample_States() START" << endl;
	 vector<int> curMarkstates,toStateVec;//fromStateVec,tempVec,
	 // vector<vector<int>> perms_frmState;
	  int bestStateInd_loc,bestStateInd;
	  vector<Double> curMarkVals_norm,curMarkVals;
	  Double dCurVal,drandomvalue,dsum_values=0.0,dinit_sum=0.0;//dfwdProb,dclstGV,dTemp,
	  int index_states_erase =0,chosen_index,cur_index;      
	  vector<Double> fwd_probsCurmark;
	
	 
	   curMarkstates = curStateSpace[n_markers-1];
	  //normalize the values so that they add upto 1	  
	  dsum_values = std::accumulate(Clst_givenGenoV[n_markers-1].begin(), Clst_givenGenoV[n_markers-1].end(),dinit_sum);
	  for(auto &kv:Clst_givenGenoV[n_markers-1])
	  {
	      curMarkVals_norm.emplace_back(kv/dsum_values);
	  }
	
	  //first sample ziM ~ p(ziM|g,v)~p(g,ZiM,v) 
	  //drandomvalue = unif(rng);
	  drandomvalue =  randZeroToOne();
	  cur_index =0;
	  drandomvalue -= curMarkVals_norm[0];
	 while((drandomvalue>0) && (cur_index<num_states-1))
	  {
		++cur_index;
		drandomvalue -= curMarkVals_norm[cur_index];		
	  }
	  curMarkVals.clear();
	 curMarkVals_norm.clear();
	  bestStateInd = curMarkstates[cur_index];
	 //cout << endl << "bestStateInd " << bestStateInd << endl;
	 chosenStates.emplace_back(bestStateInd); 

	//choose zim recursively from markers M-1,...,1
	for(int markCount=n_markers-2; markCount>=0; --markCount)
	{
		curMarkstates = curStateSpace[markCount];		
		curMarkVals=(vector<Double>(num_states));
		
		//we have value for previous marker in bestStateInd
		toStateVec = totStateSpace[bestStateInd].values;
		fwd_probsCurmark= Fwd_probs[markCount];
		
		//compute the value at this marker for all the states
		//for(int stateCount=0;stateCount< num_states;++stateCount)
		 parallel_for(int(0), num_states, [this,&fwd_probsCurmark,&toStateVec,&markCount,
			      &curMarkVals,&curMarkstates,&totStateSpace,&HmmObj](int stateCount)throw()    
		{
			Double dfwdProb,dclstGV,dTemp;
			dfwdProb =fwd_probsCurmark[stateCount];
			vector<int> fromStateVec,tempVec;
		        vector<vector<int>> perms_frmState;
			fromStateVec = totStateSpace[curMarkstates[stateCount]].values;			        
			vector_permutation(fromStateVec,tempVec,perms_frmState);		  
			dclstGV= HmmObj.compute_trans_prob_bet_clst_tuples(markCount,perms_frmState,toStateVec);	    
			dTemp = dclstGV*dfwdProb;	    
			curMarkVals[stateCount]=dTemp; 
		});
		 dsum_values = std::accumulate(curMarkVals.begin(), curMarkVals.end(),dinit_sum);
		 //cout << "dsum_values " << dsum_values << endl <<   "normalized vals" << endl;
		  for(auto &kv:curMarkVals)
		  {
		        //cout <<kv/dsum_values << " " ;
			curMarkVals_norm.emplace_back(kv/dsum_values);
		  }
		  
		  cur_index =0;
		  //drandomvalue = unif(rng);
		  drandomvalue =  randZeroToOne();
		  drandomvalue -= curMarkVals_norm[cur_index];
		  while((drandomvalue>0) && (cur_index<num_states-1))
		  {
			  ++cur_index;
			  drandomvalue -= curMarkVals_norm[cur_index];			  
		  }
		  bestStateInd = curMarkstates[cur_index];
		  //cout <<endl << markCount << " " << bestStateInd << endl;
		 chosenStates.insert(chosenStates.begin(),bestStateInd);    
		 curMarkVals.clear();
		 curMarkVals_norm.clear();
	}

  
}

void Phase::resolve_StateOrder( Chaplotypes &HmmObj,const vector<states> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates)
{
     //cout << "CHK: Phase::resolve_StateOrder() START" << endl;	
      map<pair<int,int>,Double> markerData;  
      vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
      vector<vector<int>> tostatePerms;
      vector<Double> transValues,transValues_norm;
      std ::pair<int,int> correctPair;
      int position;
      Double tempTrans,dsum_values,dinit_sum=0.0,drandomvalue;    
      vector<int> frmStateVec,toStateVec,tempVec,tempResfrm,tempResto,temp_ind;      
      
     //retain the existing order for the first marker and so start from the second marker.   
      orderedStates.emplace_back(totStateSpace[chosenStates[0]].values);
     
      //from second marker to the last marker, recursively determine the order
      for(int markCount=1;markCount< n_markers;++markCount)
      {
		vector<int> curMarkstates;
		frmStateVec = orderedStates[(markCount-1)];
		toStateVec =  totStateSpace[chosenStates[markCount]].values;      
		markerData = HmmObj.m_trans_prob_bet_clst[(markCount-1)];	    
		vector_permutation(toStateVec,tempVec,tostatePerms);
		  
		for(auto &kvp:tostatePerms)
		{
			  tempTrans = 1.0;
			  for(int kv=0;kv<n_ploidy;++kv)
			  {
			      correctPair.first = frmStateVec[kv];//kv.second;
			      correctPair.second = kvp[kv];//kv.first;
			      tempTrans *= markerData[correctPair];
			    }
			  transValues.emplace_back(tempTrans);
		    }
		//PERFORM SAMPLING
	    //normalize the values so that they add upto 1
	      dinit_sum=0.0;
	      dsum_values = std::accumulate(transValues.begin(), transValues.end(),dinit_sum);
	      for(auto &kv:transValues)
	      {
		  transValues_norm.emplace_back(kv/dsum_values);
	      }
	      position =0;
	      //drandomvalue = unif(rng);
	      drandomvalue =  randZeroToOne();
	      drandomvalue -= transValues_norm[position];
	      while(drandomvalue>0 && position < (num_states-1))
	      {
		    ++position;
		    drandomvalue -= transValues_norm[position];		 
	      }
	      //our chosen sample is at position. Now retirve ordered states
	    //  cout << endl<< "chosen sample ";
	      for(auto &kv:tostatePerms[position])
	      {
		      curMarkstates.emplace_back(kv);
		}
	    
		orderedStates.emplace_back(curMarkstates);
		Map_Valid_stateTuples.clear();
		transValues.clear();
		tostatePerms.clear();
		transValues_norm.clear();
		temp_ind.clear();
		tempResfrm.clear();
		tempResto.clear();
      }
  
  
}
  
void Phase::do_phasing(int indCount, Chaplotypes &HmmObj, 
vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase)
{
      //cout << "CHK: Phase::do_phasing() START" << endl;
      vector<vector<int>> genotype_Ind,permuted_Chaplotypes;  
      vector<int> geno_MarkerVec,cur_state,tempVec,tempResfrm,tempResto,temp_ind;
      vector< vector<std::pair<int,int>>> Map_Haplos_states;
      vector<Double> probValues,probValues_norm;
      vector<double>thetaCurMarker;
      Double perm_values,dTemphaplo,dThetaOne,dsum_values,dinit_sum=0.0,drandomvalue;    
      vector<std::pair<int,int>> tempState;
      genotype_Ind= HmmObj.m_input.m_genotypes[indCount];
      int ihaplo,position;
 
      for(int markCount=0;markCount< n_markers;++markCount)
      {
	    vector<int> phased_MarkerVec;
	    geno_MarkerVec = genotype_Ind[markCount];
	    thetaCurMarker = HmmObj.m_theta[markCount];
	    vector_permutation(geno_MarkerVec,tempVec,permuted_Chaplotypes);
	    cur_state = orderedStates[markCount];       
	    // map the state and haplo combi in unique way, for eg haplo={0, 0, 1, 1} state={2,3,3,4} map= {(0, 2),(0, 3),(1, 3),(1, 4)}
		    
	    for(auto &kvp:permuted_Chaplotypes)
	    {
		  perm_values=1.0;
		  //for(auto clst_it:kvp)
		  for(int clst_it=0;clst_it<n_ploidy;++clst_it)
		  {
			dThetaOne = thetaCurMarker[cur_state[clst_it]];
			ihaplo =   kvp[clst_it];      
			// Equation (1) in Scheet & Stephans 2006  article
			dTemphaplo = pow(dThetaOne,ihaplo);
			dThetaOne = 1.0-dThetaOne;
			ihaplo=1.0-ihaplo;  
			dTemphaplo *= pow(dThetaOne,ihaplo);
			perm_values = perm_values*dTemphaplo; 
		  }
		  //collect values from all possible permutations  
		   probValues.emplace_back(perm_values);	   
	    }   
	     
	     //PERFORM SAMPLING
	     //normalize the values so that they add upto 1
	    dsum_values = std::accumulate(probValues.begin(), probValues.end(),dinit_sum);
	    for(auto &kv:probValues)
	    {
		probValues_norm.emplace_back(kv/dsum_values);
	    }
	    position =0;
	    //drandomvalue = unif(rng);
	    drandomvalue = randZeroToOne();
	    drandomvalue -= probValues_norm[position];
	    while(drandomvalue>0 && (position<(num_states-1)))
	    {	
		  ++position;
		  drandomvalue -= probValues_norm[position];		  
	    }
	    
	    //our chosen sample is at position. Now retrieve ordered states
	   for(auto &kv:permuted_Chaplotypes[position])
	   {
		  phased_MarkerVec.emplace_back(kv);
		
	    }
	
	    
	    final_Phase.emplace_back(phased_MarkerVec);
	    permuted_Chaplotypes.clear();
	    probValues.clear();  
	    probValues_norm.clear();
	    temp_ind.clear();
	    tempResfrm.clear();
	    tempResto.clear();
      }
  // cout << "CHK: Phase::do_phasing() END" << endl;
 
}