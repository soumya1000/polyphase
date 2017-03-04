#include "samplingtf.h"

CPhasetf::CPhasetf()
{
  
}
CPhasetf::~CPhasetf()
{
  
}

void CPhasetf::resolvePhase(const  vector<vector<vector<Double>>> &Fwd_probs,vector<vector<vector<Double>>> &Clust_givenV,vector<vector<vector<int>>> &transitions,
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase,bool logging)
{
    //  cout << "CHK: CPhasetf::resolvePhase() START" << endl;
      //vector<vector<Double>> Fwd_probs_Ind,Clst_givenGenoV_Ind;
      //vector<int> chosenStates; 
      //string logFile;
      ofstream log_SFile;
      if(logging)
      {
	log_SFile.open("dip_Sampling_log.txt");
	
      }
      try
      {
	    for(int indCount =0;indCount<n_individuals;++indCount)//
	    {
	        orderedStates.emplace_back(vector<vector<int>>());
		final_Phase.emplace_back(vector<vector<int>>());
	    }
	    parallel_for(int(0),n_individuals,[this,&Fwd_probs,&Clust_givenV,&transitions,
		    &HmmObj, &totStateSpace,&orderedStates, &final_Phase](int indCount)throw()
	    {
		  vector<vector<int>> stateSpace_Ind;
		  vector<vector<int>> phased_statesInd;
		    vector<vector<Double>> Fwd_probs_Ind,Clst_givenGenoV_Ind;
                      vector<int> chosenStates; 
		    
		  Fwd_probs_Ind = Fwd_probs[indCount] ;
		  sample_States( Fwd_probs_Ind,Clust_givenV,transitions,HmmObj,totStateSpace,chosenStates);  
		  resolve_StateOrder(HmmObj,totStateSpace,chosenStates,stateSpace_Ind);  		 
		  orderedStates[indCount] =stateSpace_Ind;
		  do_phasing(indCount,HmmObj,stateSpace_Ind,phased_statesInd); 
		  final_Phase[indCount] =phased_statesInd;
		  
		/*  if(logging)
		  {
		      log_SFile << endl<<"Ind " << indCount << endl;
		      
		      log_SFile << "sample_States " <<endl;		      
		      for(auto &it:chosenStates)
			log_SFile << it << " ";
		      
		      log_SFile << endl <<"resolve_StateOrder " << endl ;		      
		      for(auto &it:stateSpace_Ind)
		      {
			  for(auto &jt:it)
			  log_SFile << jt << " ";
			 log_SFile << endl;
		      }
		      
		      log_SFile << "do_phasing " <<endl;		      
		      for(auto &it:phased_statesInd)
		      {
			  for(auto &jt:it)
			    log_SFile << jt << " ";
			  log_SFile << endl;
		      }
		  }*/
		  chosenStates.clear();
	     });
      }
      catch(const std::exception& e)
      {
	    std::cout << e.what() << '\n';
      }
   //   cout << "CHK: CPhasetf::resolvePhase() END" << endl;
}

void CPhasetf::sample_States(vector<vector<Double>> &Fwd_probs, vector<vector<vector<Double>>> &Clust_givenV,vector<vector<vector<int>>> &transitions,
		    Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,vector<int> &chosenStates)
{
        //  cout << "CHK: CPhasetf::sample_States() START" << endl;
	 vector<int> curMarkstates,toStateVec;//fromStateVec,tempVec,
	 // vector<vector<int>> perms_frmState;
	  int bestStateInd_loc,bestStateInd,state;
	  vector<Double> curMarkVals_norm,curMarkVals;
	  Double dCurVal,drandomvalue,dsum_values=0.0,dinit_sum=0.0;//dfwdProb,dclstGV,dTemp,
	  int index_states_erase =0,chosen_index,cur_index,max_index;      
	  vector<Double> fwd_probsCurmark;
	  vector<vector<Double>> Clust_givenVCurMark;
	  vector<int>::iterator stateit;
	 //normalize the values so that they add upto 1	  
	  dsum_values = std::accumulate(Fwd_probs[n_markers-1].begin(), Fwd_probs[n_markers-1].end(),dinit_sum);
	
	  for(auto &kv:Fwd_probs[n_markers-1])
	  {
	      curMarkVals_norm.emplace_back(kv/dsum_values);
	  }
	
	  //first sample ziM ~ p(ziM|g,v)~p(g,ZiM,v) 
	 
	  drandomvalue =  randZeroToOne();
	  cur_index =0;
	  drandomvalue -= curMarkVals_norm[0];
	 while((drandomvalue>0) && (cur_index<num_states-1))
	  {
		++cur_index;
		drandomvalue -= curMarkVals_norm[cur_index];		
	  }
	  bestStateInd = cur_index;
	// cout << endl << "bestStateInd " << bestStateInd << " "<< totStateSpace[bestStateInd][0]<< totStateSpace[bestStateInd][1] <<endl;
	
	 chosenStates.emplace_back(bestStateInd); 
	
	//choose zim recursively from markers M-1,...,1
	for(int markCount=n_markers-2; markCount>=0; --markCount)
	{
		//cout << endl<< "markCount " << markCount << endl;
		curMarkVals.clear();
		curMarkVals=(vector<Double>(num_states));
		curMarkVals_norm.clear();
		//we have value for previous marker in bestStateInd
		toStateVec = totStateSpace[bestStateInd];
		fwd_probsCurmark= Fwd_probs[markCount];
		Clust_givenVCurMark= Clust_givenV[markCount];
		 for(int stateCount=0;stateCount<num_states;++stateCount)
		 {
		       stateit = std::find( transitions[markCount][stateCount].begin(),transitions[markCount][stateCount].end(),bestStateInd);
		       
		       if(stateit != transitions[markCount][stateCount].end())
		      {
			    state = std::distance(transitions[markCount][stateCount].begin(),stateit);
			    curMarkVals[stateCount]= fwd_probsCurmark[stateCount]*Clust_givenVCurMark[stateCount][state];
		      }
		      else
		      {
			  curMarkVals[stateCount]=0.0;			  
		      }
		   
		}
		
		 dsum_values = std::accumulate(curMarkVals.begin(), curMarkVals.end(),dinit_sum);
		
		  for(auto &kv:curMarkVals)
		  {
		      
			curMarkVals_norm.emplace_back(kv/dsum_values);
		  }
		  
		  cur_index =0;
		  //max_index = curMarkVals_norm.size();
		  drandomvalue =  randZeroToOne();
		  drandomvalue -= curMarkVals_norm[cur_index];
		  while((drandomvalue>0) && (cur_index<num_states-1))
		  {
		      ++cur_index;
		      drandomvalue -= curMarkVals_norm[cur_index];			  
		  }
		  bestStateInd = cur_index;
		// cout <<endl <<  "bestStateInd " << bestStateInd <<" "<< totStateSpace[bestStateInd][0]<< totStateSpace[bestStateInd][1] <<endl;
		 chosenStates.insert(chosenStates.begin(),bestStateInd);    
	}
	
//	cout << "CHK: CPhasetf::sample_States() END" << endl;
}

void CPhasetf::resolve_StateOrder( Chaplotypes &HmmObj,const vector<vector<int>> &totStateSpace,
	vector<int> &chosenStates,vector<vector<int>> &orderedStates)
{
    //cout << "CHK: CPhasetf::resolve_StateOrder() START" << endl;     
      vector<Double> probValues;
      int position,numpositions,icounter;
      Double tempTrans,dsum_values,drandomvalue;    
      vector<int> toStateVec,fromstateVec,tempVec,correctOrder;      
      map< pair<int,int>,Double> trans_fromPrevMark;
      vector<vector<int>> tostatePerms;
      vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
       
     //retain the existing order for the first marker and so start from the second marker.   
      orderedStates.emplace_back(totStateSpace[chosenStates[0]]);
      
      //from second marker to the last marker, recursively determine the order      
     for(int markCount=1;markCount< n_markers;++markCount)
      {
	   trans_fromPrevMark = HmmObj.m_trans_prob_bet_clst[markCount-1];
	   fromstateVec = orderedStates[markCount-1];
	   toStateVec = totStateSpace[chosenStates[markCount]];
	   dsum_values =0.0;
	   
	   vector_permutation(toStateVec,tempVec,tostatePerms);
	   get_map_vectors(tostatePerms, fromstateVec, Map_Valid_stateTuples);
	   
	   for(auto &kvp:Map_Valid_stateTuples)
          {
		tempTrans = 1.0;
		for(auto kv:kvp)
		{
		    tempTrans *= trans_fromPrevMark[std::make_pair (kv.second,kv.first)];
		}
		probValues.emplace_back(tempTrans);
		dsum_values +=  tempTrans;
          }
	   for(auto &it:probValues)
	      it /= dsum_values; 
	   
	   position =0;
	   numpositions = probValues.size();
           drandomvalue =  randZeroToOne();
	   drandomvalue -= probValues[position];
	   while((drandomvalue>0) && (position<numpositions))
	   {
	      ++position;
	      drandomvalue -= probValues[position];			  
	   }
	   icounter =0;
	   while(icounter<n_ploidy)
	   {
	     correctOrder.emplace_back(Map_Valid_stateTuples[position][icounter].first);
	     ++icounter;
	   }
	   orderedStates.emplace_back(correctOrder);
	   correctOrder.clear();
	   probValues.clear();
	   tostatePerms.clear();
	   Map_Valid_stateTuples.clear();
      }
     //cout << "CHK: CPhasetf::resolve_StateOrder() END" << endl;	
}

void CPhasetf::do_phasing(int indCount, Chaplotypes &HmmObj, 
vector<vector<int>> &orderedStates,vector<vector<int>> &final_Phase)
{
    //cout << "CHK: CPhasetf::do_phasing() START" << endl;
      vector<vector<int>> genotype_Ind,permuted_haplotypes;  
      vector< vector<std::pair<int,int>>> Map_Haplos_states;
      vector<int> geno_MarkerVec,cur_state,phased_MarkerVec,tempVec;
      vector<Double> probValues;
      vector<double>thetaCurMarker;
      Double dTemphaplo,dsum_values,dinit_sum=0.0,drandomvalue,perm_values,dThetaOne;        
      
      int position,genosum,numpositions,ihaplo,icounter;
     
      genotype_Ind= HmmObj.m_input.m_genotypes[indCount];
      for(int markCount=0;markCount< n_markers;++markCount)
      {
	    //cout << "markCount "<<markCount << endl;
	    cur_state = orderedStates[markCount];	    
	    geno_MarkerVec = genotype_Ind[markCount];	   
	    thetaCurMarker = HmmObj.m_theta[markCount];	   
	    genosum =  std::accumulate(geno_MarkerVec.begin(),geno_MarkerVec.end(),0);
	    
	     if(genosum == 0 || genosum==n_ploidy)
	    {
	       phased_MarkerVec = geno_MarkerVec;	      
	    }
	     
	    else
	    {
	          vector_permutation(geno_MarkerVec,tempVec,permuted_haplotypes);
		  get_map_vectors(permuted_haplotypes, cur_state, Map_Haplos_states);
		   dsum_values =0.0;
		  for(auto &it:Map_Haplos_states)
		  {
			perm_values=1.0;
			for (int clst_it=0;clst_it<n_ploidy;++clst_it)
			{
			      dThetaOne = thetaCurMarker[cur_state[clst_it]];
			      ihaplo =   it[clst_it].first;			
			      // Equation (1) in Scheet & Stephans 2006  article
			      dTemphaplo = pow(dThetaOne,ihaplo);
			      dThetaOne = 1.0-dThetaOne ; 
			      ihaplo = 1.0-ihaplo;
			      dTemphaplo *= pow(dThetaOne,ihaplo);		    
			      perm_values = perm_values*dTemphaplo;	
			  }
			
			probValues.emplace_back(perm_values);
			//sum of all possible permutations     
			dsum_values +=  perm_values ;
		      }
		  
		      for(auto &it:probValues)
			it /= dsum_values;
		  
		  position =0;
		  numpositions = probValues.size();
                  drandomvalue =  randZeroToOne();
	          drandomvalue -= probValues[position];
	          while((drandomvalue>0) && (position<numpositions))
	          {
	             ++position;
	             drandomvalue -= probValues[position];			  
	           }
	          icounter =0;
		  while(icounter<n_ploidy)
	          {
	              phased_MarkerVec.emplace_back(Map_Haplos_states[position][icounter].first);
	              ++icounter;
	            }
		 
		  //cout << phased_MarkerVec[0]<< " " << phased_MarkerVec[1]  <<endl;
		  
		  Map_Haplos_states.clear();
		  permuted_haplotypes.clear();
		  probValues.clear();
	    }
	    
	   final_Phase.emplace_back(phased_MarkerVec);
	   phased_MarkerVec.clear();
      }
   //cout << "CHK: CPhasetf::do_phasing() END" << endl;
}