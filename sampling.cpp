#include "sampling.h"


Phase::Phase()
{
  
}
Phase::~Phase()
{
  
}

/*void Phase::resolvePhase(const vector<vector<vector<Double>>> &Fwd_probs,const vector<vector<vector<Double>>> &Clst_givenGenoV,
			  Chaplotypes &HmmObj,const vector<states> &totStateSpace,
		         const vector<vector<vector<int>>> &curStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)*/
void Phase::resolvePhase(const  vector<vector<vector<Double>>> &Fwd_probs,const  vector<vector<vector<Double>>> &Clst_givenGenoV,
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)
{
      cout << "CHK: Phase::resolvePhase() START" << endl;
      vector<vector<Double>> Fwd_probs_Ind,Clst_givenGenoV_Ind;
      vector<int> chosenStates;
  try
  {
	
	for(int indCount =0;indCount<n_individuals;++indCount)
	  {
		  vector<vector<int>> stateSpace_Ind;
		  vector<vector<int>> phased_statesInd;
		 
		  Fwd_probs_Ind = Fwd_probs[indCount] ;
		  Clst_givenGenoV_Ind = Clst_givenGenoV[indCount];
		 
		  sample_States( Fwd_probs_Ind,Clst_givenGenoV_Ind,HmmObj,totStateSpace,chosenStates);  
		  resolve_StateOrder(HmmObj,totStateSpace,chosenStates,stateSpace_Ind);  		 
		  orderedStates.emplace_back(stateSpace_Ind);
		  do_phasing(indCount,HmmObj,stateSpace_Ind,phased_statesInd); 
		  final_Phase.emplace_back(phased_statesInd);
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
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates)
{
	// cout << "CHK: Phase::sample_States() START" << endl;
	 vector<int> curMarkstates,toStateVec;//fromStateVec,tempVec,
	 // vector<vector<int>> perms_frmState;
	  int bestStateInd_loc,bestStateInd;
	  vector<Double> curMarkVals_norm,curMarkVals;
	  Double dCurVal,drandomvalue,dsum_values=0.0,dinit_sum=0.0;//dfwdProb,dclstGV,dTemp,
	  int index_states_erase =0,chosen_index,cur_index;      
	  vector<Double> fwd_probsCurmark;
	
	 
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
	  bestStateInd = cur_index;
	 chosenStates.emplace_back(bestStateInd); 
	
	//choose zim recursively from markers M-1,...,1
	for(int markCount=n_markers-2; markCount>=0; --markCount)
	{
		curMarkVals.clear();
		curMarkVals=(vector<Double>(num_states));
		curMarkVals_norm.clear();
		//we have value for previous marker in bestStateInd
		toStateVec = totStateSpace[bestStateInd];
		fwd_probsCurmark= Fwd_probs[markCount];
		
		//compute the value at this marker for all the states
		//for(int stateCount=0;stateCount< num_states;++stateCount)
		 parallel_for(int(0), num_states, [this,&fwd_probsCurmark,&toStateVec,&markCount,&curMarkVals,&totStateSpace,&HmmObj](int stateCount)throw()    
		{
			Double dfwdProb,dclstGV,dTemp;
			dfwdProb =fwd_probsCurmark[stateCount];
			vector<int> fromStateVec,tempVec;
		        vector<vector<int>> perms_frmState;
			fromStateVec = totStateSpace[stateCount];			        
			vector_permutation(fromStateVec,tempVec,perms_frmState);		  
			dclstGV= HmmObj.compute_trans_prob_bet_clst_tuples(markCount,perms_frmState,toStateVec);	    
			dTemp = dclstGV*dfwdProb;	    
			curMarkVals[stateCount]=dTemp; 
		});
		 dsum_values = std::accumulate(curMarkVals.begin(), curMarkVals.end(),dinit_sum);
		  for(auto &kv:curMarkVals)
		  {
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
		  bestStateInd = cur_index;
		 chosenStates.insert(chosenStates.begin(),bestStateInd);    
	}
	
	//cout << "CHK: Phase::sample_States() END" << endl;
}

void Phase::resolve_StateOrder( Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates)
{
    // cout << "CHK: Phase::resolve_StateOrder() START" << endl;	
      map<pair<int,int>,Double> markerData;  
      vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
      vector<vector<int>> tostatePerms;
      vector<Double> transValues,transValues_norm;
      std ::pair<int,int> correctPair;
      int position;
      Double tempTrans,dsum_values,dinit_sum=0.0,drandomvalue;    
      vector<int> frmStateVec,toStateVec,tempVec,tempResfrm,tempResto,temp_ind;
      
      
     //retain the existing order for the first marker and so start from the second marker.   
      orderedStates.emplace_back(totStateSpace[chosenStates[0]]);
     // cout << orderedStates[0].size() << endl;
      //from second marker to the last marker, recursively determine the order
      for(int markCount=1;markCount< n_markers;++markCount)
      {
	    vector<int> curMarkstates;
	    frmStateVec = orderedStates[(markCount-1)];
	       toStateVec =  totStateSpace[chosenStates[markCount]];      
	  markerData = HmmObj.m_trans_prob_bet_clst[(markCount-1)];
	    
	    vector_permutation(toStateVec,tempVec,tostatePerms);
	    get_map_vectors(tostatePerms, frmStateVec, Map_Valid_stateTuples);
     
	 for(auto &kvp:Map_Valid_stateTuples)
	    {
		  tempTrans = 1.0;
		  for(auto kv:kvp)
		  {
		      correctPair.first = kv.second;
		      correctPair.second = kv.first;
		      tempTrans *= markerData[correctPair];
		     // tempTrans = exp(log(tempTrans) + log(markerData[correctPair]));
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
	   for(auto &kv:Map_Valid_stateTuples[position])
	   {
		  tempResfrm.emplace_back(kv.second);
		  tempResto.emplace_back(kv.first);
	    }
	   for(auto &kv:frmStateVec)
	   {
		for(int ind=0;ind<n_ploidy;++ind)
		{
		    if(kv== tempResfrm[ind] &&(std::find(temp_ind.begin(), temp_ind.end(), ind)==temp_ind.end()))
		    {
			  temp_ind.emplace_back(ind);
			  break;
		    }
		}	
	    }
	    // the order of indices is in temp_ind
	    for(auto &kv:temp_ind)
	    {
		  curMarkstates.emplace_back(tempResto[kv]);
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
  
  //    cout << "CHK: Phase::resolve_StateOrder() END" << endl;
}

void Phase::do_phasing(int indCount, Chaplotypes &HmmObj, 
vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase)
 {
        //cout << "CHK: Phase::do_phasing() START" << endl;
      vector<vector<int>> genotype_Ind,permuted_Chaplotypes;  
      vector<int> geno_MarkerVec,cur_state,tempVec,tempResfrm,tempResto,temp_ind;
      vector< vector<std::pair<int,int>>> Map_Haplos_states;
      vector<Double> probValues,thetaCurMarker,probValues_norm;
      Double perm_values,dTemphaplo,dThetaOne,dsum_values,dinit_sum=0.0,drandomvalue;    
      vector<std::pair<int,int>> tempState;
      genotype_Ind= HmmObj.m_input.m_genotypes[indCount];
      int ihaplo,position;
         //RANDOM num generation
	/*std::mt19937_64 rng;    
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();// initialize the random number generator with time-dependent seed
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);  // initialize a uniform distribution between 0 and 1 */     
	
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
		  for(auto clst_it:kvp)
		  {
			dThetaOne = thetaCurMarker[clst_it.second];
			ihaplo =   clst_it.first;      
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
	   for(auto &kv:Map_Haplos_states[position])
	   {
		  tempResfrm.emplace_back(kv.second);
		  tempResto.emplace_back(kv.first);
	    }
	   for(auto &kv:cur_state)
	   {
		for(int ind=0;ind<n_ploidy;++ind)
		{
		      if(kv== tempResfrm[ind] &&(std::find(temp_ind.begin(), temp_ind.end(), ind)==temp_ind.end()))
		      {
			    temp_ind.emplace_back(ind);
			    break;
		      }
		}	
	    }
	    // the order of indices is in temp_ind
	    for(auto &kv:temp_ind)
	    {
		  phased_MarkerVec.emplace_back(tempResto[kv]);
	    }
	    
	    final_Phase.emplace_back(phased_MarkerVec);
	    permuted_Chaplotypes.clear();
	    Map_Haplos_states.clear();
	    probValues.clear();  
	    probValues_norm.clear();
	    temp_ind.clear();
	    tempResfrm.clear();
	    tempResto.clear();
      }
  
}


void Phase::sample_States_OLD(vector<vector<Double>> &Fwd_probs, vector<vector<Double>> &Clst_givenGenoV, 
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates)
{
  vector<int> curMarkstates,fromStateVec,tempVec,toStateVec;
  vector<vector<int>> perms_frmState;
  int bestStateInd_loc,bestStateInd;
  vector<Double> curMarkValues;
  Double dfwdProb,dclstGV,dCurVal,dTemp;
  
  //first sample ziM ~ p(ziM|g,v)~p(g,ZiM,v) 
  bestStateInd_loc = distance(begin(Clst_givenGenoV[n_markers-1]), max_element(begin(Clst_givenGenoV[n_markers-1]), end(Clst_givenGenoV[n_markers-1])));  
  //bestStateInd = curStateSpace[n_markers-1][bestStateInd_loc];
  bestStateInd = bestStateInd_loc;
  chosenStates.emplace_back(bestStateInd);
  dCurVal = Clst_givenGenoV[n_markers-1][bestStateInd_loc];  
  
  //cout << n_markers-1 << " " << bestStateInd_loc << " " << bestStateInd <<endl;
  //choose zim recursively from markers M-1,...,1
  for(int markCount=n_markers-2; markCount>=0; --markCount)
  {
	  //we have value for previous marker in bestStateInd
	  fromStateVec = totStateSpace[bestStateInd];
	  //curMarkstates = curStateSpace[markCount];
	  vector_permutation(fromStateVec,tempVec,perms_frmState);        
	  //compute the value at this marker for all the states
	  for(int stateCount=0;stateCount< num_states;++stateCount)
	  {
		  dfwdProb = Fwd_probs[markCount][stateCount];
		// toStateVec = totStateSpace[curMarkstates[stateCount]].values;
		  toStateVec = totStateSpace[stateCount];		  
		  dclstGV= HmmObj.compute_trans_prob_bet_clst_tuples(markCount,perms_frmState,toStateVec);		  
		  dTemp = dCurVal*dfwdProb;
		  dTemp *= dclstGV;
		  curMarkValues.emplace_back(dTemp);  
	  }
	
	  //choose the state that has maximum value at this marker
	  bestStateInd_loc = distance(begin(curMarkValues), max_element(begin(curMarkValues), end(curMarkValues))); 
	  dCurVal *= curMarkValues[bestStateInd_loc];      
	  bestStateInd = bestStateInd_loc;	  
	  chosenStates.emplace(chosenStates.begin(),bestStateInd); 
	  curMarkValues.clear();
	  perms_frmState.clear();      
  }
     

}

void Phase::resolve_StateOrder_OLD( Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates,vector<vector<int>> &orderedStates)
{
  //get_best_orderToState(int frmMarker,vector<int>&frmState,vector<vector<int>>&tostatePerms,vector<int>&tobestState)
  //  cout << "CHK: Phase::resolve_StateOrder() START" << endl;
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
   orderedStates.emplace_back(totStateSpace[chosenStates[0]]);
  
   //from second marker to the last marker, recursively determine the order
   for(int markCount=1;markCount<n_markers;++markCount)
   {
      vector<int> curMarkstates;
      
      frmStateVec = orderedStates[markCount-1];
      //toStateVec =  totStateSpace[chosenStates[markCount]].values; 
      toStateVec =  totStateSpace[chosenStates[markCount]]; 
     
      markerData = HmmObj.m_trans_prob_bet_clst[markCount-1];
      vector_permutation(toStateVec,tempVec,tostatePerms);
      get_map_vectors(tostatePerms, frmStateVec, Map_Valid_stateTuples);
     
      for(auto &kvp:Map_Valid_stateTuples)
      {
	    tempTrans = 1.0;
	    for(auto kv:kvp)
	    {
		correctPair.first = kv.second;
		correctPair.second = kv.first;
		tempTrans = exp(log(tempTrans) + log(markerData[correctPair]));
	    }
	    transValues.emplace_back(tempTrans);
	 }
   
      position =distance(begin(transValues), max_element(begin(transValues), end(transValues)));
      tempState = Map_Valid_stateTuples[position];
 
      for(auto &kv:tempState)
      {
	  curMarkstates.emplace_back(kv.first);
      }
     
      orderedStates.emplace_back(curMarkstates);
      Map_Valid_stateTuples.clear();
      transValues.clear();
      tostatePerms.clear();
   }
  
}

void Phase::do_phasing_OLD(int indCount, Chaplotypes &HmmObj, 
vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase)
{
   //cout << "CHK: Phase::do_phasing() START" << endl;
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
 //for(auto &kvp:final_Phase)
 //  {
     //for(auto kv:kvp)
      //cout << kv << ",";     
     // cout << "  ";
   //}
 //  cout << endl;
  //  cout << "CHK: Phase::do_phasing() END" << endl;
}