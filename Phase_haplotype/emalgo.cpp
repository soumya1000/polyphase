
#include "emalgo.h"


//*************************************************************//
// CLASS EM_HMM_GENOTYPE
//*************************************************************//

em_hmm_genotype::em_hmm_genotype(string ipFileName)
{
      cout << "em_hmm_genotype::em_hmm_genotype()" << endl;
      m_params_file.open ("Phase_param.txt"); 
      //m_GSL_file.open("GSL_param.txt"); 
      initialise_model_param(ipFileName);
      get_total_stateSpace();   
      get_current_stateSpace(); 
     // update_HMM_param();
      //log_param(); 
 }

em_hmm_genotype::~em_hmm_genotype()
{
   cout << "em_hmm_genotype::~em_hmm_genotype()" << endl;
   m_params_file.close();   
   //m_GSL_file.close();
}

void em_hmm_genotype::initialise_model_param(string ipFileName)
{
  m_func_Scaler =0; // this should be zero when ever the param-Opt is started.
  m_hap_data.initialise_param(ipFileName); 
 // m_hap_data.log_param_hap();  

}
void em_hmm_genotype::get_total_stateSpace()
{
     //cout << "em_hmm_genotype::get_total_stateSpace() START" << endl;
      //gen_combi_states(n_ploidy,n_clusters,m_total_states); 
      gen_combi_states(m_total_states); 
      m_stateweights= std::vector<Double>(m_total_states.size());
    std::transform( m_total_states.begin(), m_total_states.end(), m_stateweights.begin(),
        [] (states const& ms){ return ms.weight;});
  
  // cout << "em_hmm_genotype::get_total_stateSpace() END" << endl;
}

void em_hmm_genotype::update_HMM_param()
{
   //cout<< "CHK em_hmm_genotype::update_HMM_param() START" << endl; 
     clear_variables();
         
      m_hap_data.compute_trans_prob_bet_clst();
      compute_jump_prob();        

      for(int indCount=0;indCount<n_individuals;++indCount)
	{
	      m_clst_given_v .emplace_back(vector<vector<Double>>());
	      m_geno_given_clst_v.emplace_back(vector<vector<Double>>());   
	      m_clst_given_geno_v.emplace_back(vector<vector<Double>>());
	      m_fwd_probs.emplace_back(vector<vector<Double>> ());
	      m_bckwd_probs.emplace_back(vector<vector<Double>> ());      
	}
      
	for(int indCount=0;indCount<n_individuals;++indCount)
	{
	      for(int marker=0;marker<n_markers;++marker)
	      {
		   m_fwd_probs[indCount].emplace_back(vector<Double>(num_states));
		  m_bckwd_probs[indCount].emplace_back(vector<Double>(num_states));
	      }
	}
      
try
      {
	   compute_clust_given_v();  	
	    compute_geno_given_clst_v();	
	    compute_fwd_prob();  	
	    compute_bckwd_prob(); 
	    compute_clst_given_geno_v(); 
	    
      }
      catch(const std::exception& e)
      {
	    std::cout << e.what() << '\n';
      }
      
   //cout<< "CHK em_hmm_genotype::update_HMM_param() END" << endl;
}

//p(Jim=j|r)
void em_hmm_genotype::compute_jump_prob()
{
  // cout<< "CHK em_hmm_genotype::compute_jump_prob() START" << endl;
   
   Double jumpProbTemp,jumpProb,dTemp;
   int iscaledTimes;
   Double dCopyRecomb; 
   
   for(int jCount=0; jCount < n_markers-1; jCount++)
   {
	vector<Double> markerdata; 
	
	dTemp= m_hap_data.m_recombinations[jCount]*m_hap_data.m_input.m_physical_distances[jCount];
	  
	jumpProbTemp = exp(-1*dTemp);
	
	
	for(int iCount=0;iCount <= n_ploidy;iCount++)
	{
	    jumpProb = binomialProb(n_ploidy,iCount,jumpProbTemp);
	    
	    markerdata.push_back(jumpProb); 
	}
        
        m_jump_prob[jCount] = markerdata;  
    } 
  
 // cout << "em_hmm_genotype::compute_jump_prob() END" << endl;
}

void em_hmm_genotype::get_current_stateSpace()
{
   //cout << "em_hmm_genotype::get_current_stateSpace() START" << endl;
      //we sample separate state space for different individuals
      for(int indCnt=0;indCnt<n_individuals;++indCnt)
      {
		vector<vector<int>> indData;
		//cout << endl << indCnt<< endl;
		for(int markCount=0; markCount < n_markers; markCount++)   
		{
			//cout << endl<<markCount << endl;
			vector<int> markerData;	      
			//sample indices of 10 states according to theirs weights randomly
			sample(m_stateweights,markerData,num_states);	
			//for(auto kt : markerData)
			  //cout << kt << " " ;
			indData.emplace_back(markerData);
		}	    
		m_current_states.emplace_back(indData);
	}
  //cout << "em_hmm_genotype::get_current_stateSpace() END" << endl;
}
void em_hmm_genotype::compute_clust_given_v() 
{
    //cout<< "CHK em_hmm_genotype::compute_clust_given_v() START" << endl;
  
      try
     {
	    //for initial marker 0 
	      vector<int> cur_states, prev_states;
	      for(int indCount=0; indCount < n_individuals; indCount++)   
	      {	 
			vector<Double>initialstateprobs(num_states);
			cur_states = m_current_states[indCount][0];
			parallel_for(int(0), num_states, [this,&initialstateprobs,&cur_states](int stateCount)throw()
			{
				vector<int> toStVec,tempVec;
				vector<vector<int>> fromPerms;
				Double dFinalValue =1.0;
				toStVec = m_total_states[cur_states[stateCount]].values;
				vector_permutation(toStVec,tempVec,fromPerms);
				  
				for(auto  &iter:toStVec)
				{
				      dFinalValue = dFinalValue * m_hap_data.m_alpha[0][iter];
				}
				dFinalValue = dFinalValue * fromPerms.size();
				initialstateprobs[stateCount]=(dFinalValue);	 
				//fromPerms.clear();
			  });
			  m_clst_given_v[indCount].emplace_back(initialstateprobs);
		
			//for markers 1 to n_markers  
			for(int markCnt=1; markCnt < n_markers;++markCnt)
			{
				vector<Double> markerData(num_states);
				cur_states = m_current_states[indCount][markCnt];
				prev_states = m_current_states[indCount][markCnt-1];
				parallel_for(int(0), num_states, [this,&markerData,&markCnt,&prev_states,&cur_states](int stateCount)throw()
				{
					vector<int> toStVec;			    
					toStVec  = m_total_states[cur_states[stateCount]].values;
					
					combinable<Double> dFinalValue([]() { return 0.0; }); 
					// for(int fromCount=0;fromCount< num_states;++fromCount)//num_states
					  parallel_for(int(0), num_states, [this,&toStVec,&markCnt,&dFinalValue,&prev_states](int fromCount)throw()
					  {
						  vector<int> frmStVec,tempVec;
						  vector<vector<int>> fromPerms;
						  Double dTemp;
						  frmStVec = m_total_states[prev_states[fromCount]].values;
						  //get the permutations of the from state vector
						  vector_permutation(frmStVec,tempVec,fromPerms);
						  dTemp =  m_hap_data.compute_trans_prob_bet_clst_tuples(markCnt-1,fromPerms,toStVec);
						  dFinalValue.local() += dTemp;
						  fromPerms.clear();		      
					  });					  				
					  markerData[stateCount]=dFinalValue.combine(plus<Double>());
					  dFinalValue.clear();
				  });
				// m_clst_given_v[markCnt]=markerData;   
				  m_clst_given_v[indCount].emplace_back(markerData);
			}
	      }
      
 }
  catch(const std::exception& e)
  {
	std::cout << e.what() << '\n';
  }
   //   cout<< "CHK em_hmm_genotype::compute_clust_given_v() END" << endl
}
//p(zim|v)


void em_hmm_genotype::compute_geno_given_clst_v()
{
  //cout<< "CHK em_hmm_genotype::compute_geno_given_clst_v() START" << endl;
 
  vector<int> currentMark_states,temp;
  vector<double> thetaCurMarker;
  vector<Double>  clustCurMarker;
  vector<int> haplo_vector;
  vector<vector<int>> genoCurInd, permuted_haplotypes,indStates;

  vector< vector<std::pair<int,int>>> Map_Haplos_states;
 try
  {
         for(int indCount=0;indCount<n_individuals;++indCount)
	 {
	      genoCurInd = m_hap_data.m_input.m_genotypes[indCount];   
	
	      for(int markCount=0; markCount < n_markers; ++markCount)  
	      { 
		    vector<Double> markerData(num_states);
		    thetaCurMarker = m_hap_data.m_theta[markCount];
		      clustCurMarker = m_clst_given_v[indCount][markCount];
		      currentMark_states = m_current_states[indCount][markCount];
		    //get the haplo vector.For eg if geno is 0 for tetraploids,  haplo_vector= {0, 0, 0, 0} as given in the inout file
		    haplo_vector = genoCurInd[markCount];  
		    vector_permutation(haplo_vector,temp,permuted_haplotypes);
		   
		    parallel_for(int(0), num_states, [this,&permuted_haplotypes,&markerData,&thetaCurMarker,&clustCurMarker,&currentMark_states](int stateCount)throw()
		    {
			  long int ihaplo;
			  Double final_value =0.0,dTemphaplo,perm_values,dThetaOne; 
			  vector<int>  cur_state;
			  vector< vector<std::pair<int,int>>> Map_Haplos_states;
			 
			  cur_state = m_total_states[currentMark_states[stateCount]].values;
		
			  for(auto &kvp:permuted_haplotypes)
			  {
			      perm_values=1.0;
			    
				for (int clst_it=0;clst_it<n_ploidy;++clst_it)
				{
				      dThetaOne = thetaCurMarker[cur_state[clst_it]];
				      ihaplo =   kvp[clst_it];			      
				      // Equation (1) in Scheet & Stephans 2006  article
				      dTemphaplo = pow(dThetaOne,ihaplo);
				      dThetaOne = 1.0-dThetaOne ; 
				      ihaplo = 1.0-ihaplo;
				      dTemphaplo *= pow(dThetaOne,ihaplo);		    
				      perm_values = perm_values*dTemphaplo;		 
			      }
			    //sum of all possible permutations     
			    final_value +=  perm_values ;	    
			}
		       final_value *= clustCurMarker[stateCount];
			markerData[stateCount]=(final_value);
		    });
		    permuted_haplotypes.clear();  
		    m_geno_given_clst_v[indCount].emplace_back(markerData);
		}
	 }
      
  }
  catch(const std::exception&e)
  {
     std::cout << e.what() << '\n';
  }
    //cout << "em_hmm_genotype::compute_geno_given_clst_v() END" << endl;

}
 void  em_hmm_genotype::util_hmm_fwd(int individual, int marker,int state,int njumps,vector<vector<int>> &jumpData, 
			  vector<Double>&alpha_data)
{
    // cout<< "CHK em_hmm_genotype::util_hmm_fwd() START" << endl;
      vector<double> AlphaMarker;
      int state_jump;
      Double dalpha_temp;
      vector<int> input_comb,copy_Curstate;
      
      vector<int> prev_markerStates = m_current_states[individual][marker-1];
      
      // WE FIRST work to obtain the jumpData
      vector<int> curState= m_total_states[ m_current_states[individual][marker][state]].values;
      vector<int> gotTemp;  
      vector< vector<int>> T,Tfinal;
      vector< vector<int>> ::iterator itTemp;
      // fetch all the possible combinations of clusters with given number of jumps in current state
      vector_combination(curState,gotTemp, 0, njumps,T);

      for(int i=0;i<n_clusters;++i)
      {
	    input_comb.emplace_back(i);
      }
      vector<int>::iterator pos;   
      std::vector<vector<int>> state_values(num_states);
	std::transform( prev_markerStates.begin(), prev_markerStates.end(), state_values.begin(),
	    [this] (int const& ms){ return m_total_states[ms].values;});
      for(auto &it:T)
      {
	      copy_Curstate = curState;
	      vector<int> jumpDatacur;
	      
	      //retain the clusters within the state, where jump does not occur
		for(auto &jt:it)
		{
			pos=std::find(copy_Curstate.begin(), copy_Curstate.end(), jt);
			if(pos!= copy_Curstate.end())
			      copy_Curstate.erase(pos);
		}
		// fetch all possible clusters in positions within the state, where jump occurs
	      vector_combination(input_comb,gotTemp, 0, njumps,Tfinal);
		
		// combine both the above values to retieve the required state-Ids
		for(auto &iter:Tfinal)
		{
			for(auto &jt:copy_Curstate)
			    iter.emplace_back(jt);
			
			    std::sort(iter.begin(),iter.end());
			itTemp = std::find(state_values.begin(),state_values.end(), iter);
		      
			if(itTemp!=state_values.end())
			{
			    state_jump = std::distance( state_values.begin(),itTemp);
			    jumpDatacur.emplace_back(state_jump);			
			}
		}
		jumpData.emplace_back(jumpDatacur);
		Tfinal.clear();
      }
	    
    // Now we compute the alpha Values accordingly;
    
	    AlphaMarker =  m_hap_data.m_alpha[marker];	
	    
	    for(auto &it:T)
	    {
		  vector<Double> alpha_cur;
		  dalpha_temp = 1.0;
		
		  for(auto jt:it)
		      dalpha_temp*=AlphaMarker[jt];
		
		  alpha_data.emplace_back(dalpha_temp);
	      
	    }   
      // cout<< "CHK em_hmm_genotype::util_hmm_fwd() END" << endl;
} 

void em_hmm_genotype::util_hmm_bckwd(int individual,int marker,int state,int njumps,vector<vector<int>> &jumpData,
			vector<vector<Double>>&alpha_data)
{
      vector<double> AlphaMarker;
      int state_jump;
      Double dalpha_temp;
      vector<int> input_comb,copy_Curstate;
      vector< vector<int>> ::iterator itTemp;
      
      vector<int> next_markerStates = m_current_states[individual][marker+1];
      // WE FIRST work to obtain the jumpData
      vector<int> curState= m_total_states[ m_current_states[individual][marker][state]].values;
      vector<int> gotTemp;  
      vector< vector<int>> T,Tfinal;
      vector<vector< vector<int>>> Talpha;
	// fetch all the possible combinations of clusters with given number of jumps in current state
      vector_combination(curState,gotTemp, 0, njumps,T);

      
      for(int i=0;i<n_clusters;++i)
      {
	  input_comb.emplace_back(i);
      }
      vector<int>::iterator pos;
      
      std::vector<vector<int>> state_values(num_states);
      std::transform( next_markerStates.begin(), next_markerStates.end(), state_values.begin(),
	    [this] (int const& ms){ return m_total_states[ms].values;});
	
      for(auto &it:T)
      {
	      copy_Curstate = curState;
	      vector<int> jumpDatacur;
	      vector< vector<int>> Talpha_temp;
	      //retain the clusters within the state, where jump does not occur
		for(auto &jt:it)
		{
		      pos=std::find(copy_Curstate.begin(), copy_Curstate.end(), jt);
		      if(pos!= copy_Curstate.end())
			  copy_Curstate.erase(pos);
		}
		
		// fetch all possible clusters in positions within the state, where jump occurs
		vector_combination(input_comb,gotTemp, 0, njumps,Tfinal);
		
		// combine both the above values to retieve the required state-Ids
		for(auto &iter:Tfinal)
		{
			Talpha_temp.emplace_back(iter);
			for(auto &jt:copy_Curstate)
			    iter.emplace_back(jt);
			
			std::sort(iter.begin(),iter.end());
			itTemp = std::find( state_values.begin(), state_values.end(), iter);
			
			if(itTemp!=state_values.end())
			{
				state_jump = std::distance( state_values.begin(),itTemp);
				jumpDatacur.emplace_back(state_jump);			
			}	      
		}
		jumpData.emplace_back(jumpDatacur);
		Talpha.emplace_back(Talpha_temp);
		Tfinal.clear();
      }
	  
 // Now we compute the alpha Values accordingly;
 
	AlphaMarker =  m_hap_data.m_alpha[marker+1];	
	
	 for(auto &it:Talpha)
	{
	    vector<Double> alpha_cur;	    	
	    
	    for(auto jt:it)
	    {
		dalpha_temp = 1.0;
		for(auto &kt:jt)
		{
		    dalpha_temp*=AlphaMarker[kt];
		}
		 alpha_cur.emplace_back(dalpha_temp);
	    }
	     alpha_data.emplace_back(alpha_cur);	  
	}    
  //cout<< "CHK em_hmm_genotype::util_hmm_bckwd() END" << endl;
}

void em_hmm_genotype::compute_fwd_prob()
{
	//cout<< "CHK em_hmm_genotype::compute_fwd_prob()" << endl;  
      vector<Double> genoGivenClst_cur; 
      vector<vector<Double>>  genoCurInd;
      Double dTempFwdProb; 
      vector<double> AlphaCurMarker;
      vector<int> cur_markerStates,prevMark_states;
      vector<Double> jumpProbPrevMarker;
      vector<vector<Double>> FwdProbCurInd;
      
      AlphaCurMarker = m_hap_data.m_alpha[0];
      for(int indCount=0;indCount<n_individuals;++indCount)
      {
	      // indStates = m_current_states[indCount];
	      genoGivenClst_cur = m_geno_given_clst_v[indCount][0];
	      cur_markerStates = m_current_states[indCount][0];
	      parallel_for(int(0), num_states, [this,&genoGivenClst_cur,&AlphaCurMarker,&indCount,&cur_markerStates](int stateCount)throw()    
	      {
		    vector<int> curState= m_total_states[cur_markerStates[stateCount]].values;
		    Double dTempFwdProb =1.0;
		    //for(vector<int>::iterator j = curState.values.begin();j != curState.values.end();j++)
		    for(vector<int>::iterator j = curState.begin();j != curState.end();++j)
		    {
			  dTempFwdProb *= AlphaCurMarker[(*j)];
		    }
		    dTempFwdProb *=  genoGivenClst_cur[stateCount];
		    m_fwd_probs[indCount][0][stateCount]=(dTempFwdProb);
		});
	      
	     for(int markCnt=1; markCnt < n_markers ;markCnt++)
	      {
		    jumpProbPrevMarker = m_jump_prob[markCnt-1];
		    AlphaCurMarker = m_hap_data.m_alpha[markCnt];
		    cur_markerStates = m_current_states[indCount][markCnt];
		    prevMark_states = m_current_states[indCount][markCnt-1];
		    
		   parallel_for(int(0), num_states, [this,&jumpProbPrevMarker,&AlphaCurMarker,
				&markCnt,&cur_markerStates,&prevMark_states,&indCount](int toStCnt)throw()    
		     {
			        
				Double jumpVal,dsumalphaVal,dTemp,dFwdProbVal = 0.0;		    
				vector<vector<vector<int>>> JumpDatavec;
				vector<vector<Double>> alphaValTostates;
				vector<vector<int>> states_jump;
				vector<Double> alpha_jump;
				
				vector<int> toStVec=  m_total_states[cur_markerStates[toStCnt]].values;
				vector<Double>FwdProbPrevMarker = m_fwd_probs[indCount][markCnt-1];	
								
				Double dalphaValstate=1.0;	  
				//required for case where njumps=n_ploidy				
				 for(int ncount=0;ncount<n_ploidy;++ncount)
				 {
					  dalphaValstate *= AlphaCurMarker[toStVec[ncount]];
				  }				
				
				//num of jumps is zero
				auto position =std::find(prevMark_states.begin(), prevMark_states.end(), toStCnt);			      
				 if(position != prevMark_states.end())
					dFwdProbVal = FwdProbPrevMarker[std::distance(prevMark_states.begin(),position)] * jumpProbPrevMarker[0];// 
				//for jumps from 1 to n_ploidy, get jumpProbability 
				//and corresponding alphaval contributionbased on num jumps
				for(int jumpCount=1; jumpCount < n_ploidy;++jumpCount)
				{
				      vector<vector<int>> curJumpData;
				      vector<Double>alphaValTostates_Temp;
				      util_hmm_fwd(indCount,markCnt,toStCnt,jumpCount,curJumpData,alphaValTostates_Temp);			   
				      JumpDatavec.emplace_back(curJumpData);	
				      alphaValTostates.emplace_back(alphaValTostates_Temp);
				}
				
				 
				 //num of jumps is zero
				 dFwdProbVal = FwdProbPrevMarker[toStCnt] * jumpProbPrevMarker[0];//
				  
				  for(int jumpCount=1; jumpCount < n_ploidy;++jumpCount)
				  {
					    jumpVal = (jumpProbPrevMarker[jumpCount])/binomialCoef(n_ploidy,jumpCount);	
					    states_jump = JumpDatavec[jumpCount-1];		       
					    alpha_jump = alphaValTostates[jumpCount-1];        
					    //get alpha values & values of fwd probs from prev marker
					    dsumalphaVal=0 ;   	      
					    for(int ncount=0;ncount<states_jump.size();++ncount)
					    {
						  dTemp=0;		
						  for(int ncomb_count=0;ncomb_count<states_jump[ncount].size();++ncomb_count)
						  {
							dTemp+= FwdProbPrevMarker[states_jump[ncount][ncomb_count]];
						  }
						  dTemp *= alpha_jump[ncount];
						  dsumalphaVal +=(dTemp);
					    }
					    dFwdProbVal += (dsumalphaVal*jumpVal);					 
				     }
				      
				      // for jump= n_ploidy
				       dTemp=0;	
				        for(int ncount=0;ncount<num_states;++ncount)
					{
					      dTemp+= FwdProbPrevMarker[ncount];
					 }
					 dFwdProbVal += (dTemp *jumpProbPrevMarker[n_ploidy]*dalphaValstate);		 
					 dFwdProbVal *= m_geno_given_clst_v[indCount][markCnt][toStCnt];
					 m_fwd_probs[indCount][markCnt][toStCnt]=dFwdProbVal;	
		     });
	      }//MARKER ENDS
	     
      }//ind ends
      
 
  
 //cout << "em_hmm_genotype::compute_fwd_prob() END" << endl;
       
}

void em_hmm_genotype::compute_bckwd_prob()
{
  
  //cout << "em_hmm_genotype::compute_bckwd_prob() START" << endl;
    
   vector<Double> BckwrdProbNxtMarker;
   vector<Double> jumpProbNxtMarker;
   vector<double> AlphaNXtMarker;
   vector<int> cur_markerStates,next_markerStates;
   //The backward prob at the Final marker is 1 
   for(int indCount=0;indCount<n_individuals;++indCount)
   {
	  for(int i=0;i < num_states;++i)
	  {
		m_bckwd_probs[indCount][n_markers-1][i]=(1.0);
	   }	    
	  // for remaining markers
	  for(int markCnt = n_markers-2 ; markCnt >= 0 ;--markCnt)
	  {
		     jumpProbNxtMarker = m_jump_prob[markCnt];
		    AlphaNXtMarker = m_hap_data.m_alpha[markCnt+1];
		    cur_markerStates = m_current_states[indCount][markCnt];
		    next_markerStates = m_current_states[indCount][markCnt+1];
		   parallel_for(int(0), num_states, [this,&jumpProbNxtMarker,&AlphaNXtMarker,&cur_markerStates,&next_markerStates,
				&markCnt,&indCount](int toStCnt)throw()    
		    {
			      vector<vector<int>> states_jumpTemp;
			      vector<vector<Double>> alpha_jumpTemp;
			      vector<int> toStVec,curJumpstates;
			      Double jumpVal,dsumalphaVal,dTemp,dBwdProbVal = 0.0;		
			      toStVec = m_total_states[cur_markerStates[toStCnt]].values;			      
			      vector<vector<vector<int>>> JumpDatavec;
			      vector<vector<vector<Double>>> alphaValTostates;
			      
			      vector<Double>BckwrdProbNxtMarker =m_bckwd_probs[indCount][markCnt+1];
			      vector<Double>geno_given_clst_v = m_geno_given_clst_v[indCount][markCnt+1];    
			      //for jumps from 1 to n_ploidy, get jumpProbability 
			      //and corresponding alphaval contributionbased on num jumps
			      for(int jumpCount=1; jumpCount < n_ploidy;++jumpCount)
			      {
				    vector<vector<int>> curJumpData;
				    vector<vector<Double>>alphaValTostates_Temp;
				    util_hmm_bckwd(indCount,markCnt,toStCnt,jumpCount,curJumpData,alphaValTostates_Temp);			   
				    JumpDatavec.emplace_back(curJumpData);	
				    alphaValTostates.emplace_back(alphaValTostates_Temp);
			      }
			
			      //num of jumps is zero
			      auto position =std::find(next_markerStates.begin(), next_markerStates.end(), toStCnt);
			      if(position != next_markerStates.end())
			      {
				    dBwdProbVal = BckwrdProbNxtMarker[std::distance(std::begin(next_markerStates),position)]
				    *jumpProbNxtMarker[0]* geno_given_clst_v[toStCnt];
			      }
			       
			      //for jumps 1 to n
			      for(int jumpCount=1; jumpCount < n_ploidy;++jumpCount)
			      {
					jumpVal = (jumpProbNxtMarker[jumpCount])/binomialCoef(n_ploidy,jumpCount);	
					states_jumpTemp = JumpDatavec[jumpCount-1];	
					alpha_jumpTemp= alphaValTostates[jumpCount-1];
					dsumalphaVal =0;
					
					//get alpha values & values of bwd probs from Nxt marker
					for(int ncount=0;ncount<states_jumpTemp.size();++ncount)
					{
					      dTemp=0;	
					      curJumpstates = states_jumpTemp[ncount];
					      for(int ncomb_count=0;ncomb_count<curJumpstates.size();++ncomb_count)
					      {
						    dTemp+=( BckwrdProbNxtMarker[states_jumpTemp[ncount][ncomb_count]]*
						    geno_given_clst_v[states_jumpTemp[ncount][ncomb_count]] * alpha_jumpTemp[ncount][ncomb_count]);
					      }
					      dsumalphaVal +=   dTemp;
					}
					dBwdProbVal +=  (dsumalphaVal*jumpVal);
			      }
		  
			      //num of jumps is n_ploidy
			    dsumalphaVal =0;
			    for(int ncount=0;ncount < next_markerStates.size();++ncount)
			    {
				  dTemp = 1.0;
				  for(int icount = 0;icount<n_ploidy;++icount)
				    dTemp *=AlphaNXtMarker[icount];
				
				  dsumalphaVal += (dTemp *geno_given_clst_v[ncount]*BckwrdProbNxtMarker[ncount]);			
			    }
			    dsumalphaVal *=  jumpProbNxtMarker[n_ploidy];
			    dBwdProbVal +=  dsumalphaVal;
			    m_bckwd_probs[indCount][markCnt][toStCnt]=dBwdProbVal;	 
		    });
	  }
    }
  

 //cout << " em_hmm_genotype::compute_bckwd_prob() END " << endl;
  
}

Double em_hmm_genotype::compute_geno_given_v(int ind, int mark)
{
   // cout << "em_hmm_genotype::compute_geno_given_v() Start" << endl;
  
   Double dTempmark=0.0;
    int stateIndex;
    vector<Double > genoGivenClst_cur,fwdProb;
  
   // int iscaledTimes = 0;
    
    if(mark >=n_markers)
    {
      return(-1.0);
    }
    
    fwdProb = m_fwd_probs[ind][mark];
    
    for(std::vector<Double>::iterator i=fwdProb.begin();i!=fwdProb.end();i++) 
    {
      dTempmark = dTempmark + ((*i));
    }
   
  //  cout << "em_hmm_genotype::compute_geno_given_v()  End" << endl; 
    
    return(dTempmark); 
}

void em_hmm_genotype::compute_clst_given_geno_v()
{
 
  //cout<< "CHK compute_clst_given_geno_v Start" << endl;  
  vector<Double> FwdProbCurMark, BckwrdProbCurMark;  
  vector<int> currentMark_states;
  Double dTemp;
  
  for(int indCnt = 0; indCnt < n_individuals; ++indCnt)
  {
        
	  for(int markCnt=0; markCnt < n_markers ;++markCnt)
	  {
		  vector<Double> markerData(num_states);
		  FwdProbCurMark = m_fwd_probs[indCnt][markCnt];
		  BckwrdProbCurMark = m_bckwd_probs[indCnt][markCnt];
		  currentMark_states =  m_current_states[indCnt][markCnt];
	    
		// for(int stateCount=0; stateCount < num_states; stateCount++)
		  parallel_for(int(0), num_states, [this,&FwdProbCurMark,&BckwrdProbCurMark,&markerData,&currentMark_states](int stateCount)throw()    
		  {
			  vector<int>  stateVec = m_total_states[currentMark_states[stateCount]].values;	      
			  markerData[stateCount]= FwdProbCurMark[stateCount]*BckwrdProbCurMark[stateCount];
		  });
		  
		    m_clst_given_geno_v[indCnt].emplace_back(markerData);
	      }
   }
  

  //cout<< "CHK compute_clst_given_geno_v End" << endl;
    
}
  
void em_hmm_genotype::log_param()
{
         m_hap_data.log_param_hap(); 
	    
	  int marker =0, indcnt=0,state;
      
	  m_params_file << "..............................statespace........................................" << endl;    
	  
	  for(auto  &it:m_current_states)
	  {
	      marker =0;
	      m_params_file << "Ind: " << indcnt << " " ;
	      for(auto jt : it)
	      {   
		m_params_file << "marker: " << marker << " " ;
		  for(auto kt : jt)
		    m_params_file << kt << " " ;
		++marker;
		m_params_file << endl;
	      }
	      ++indcnt;
	      m_params_file << endl;
	      
	  }
	
	  m_params_file << "..............................JumpProb........................................" << endl;
	  for(int jCount=0; jCount < n_markers-1; jCount++)
	  {
	      for(int iCount=0;iCount <= n_ploidy;iCount++)
	      {
		m_params_file << " " << m_jump_prob[jCount][iCount];
	      }
	      m_params_file << endl;
	  }
	  
	    m_params_file << endl<<"***** compute_clust_given_v*****"  << endl;
	  indcnt=0;
	  for(auto  &it:m_clst_given_v)
	  { 
		m_params_file<<endl << "Ind: " << indcnt << endl ;
		marker=0;
		for(auto jt : it)
		{   
		      m_params_file<<endl << "marker: " << marker << endl ;
		      for(auto kt : jt)
			    m_params_file << " "  <<kt<< " " ;
		      ++marker;
		}
		++indcnt;
	  } 
	
	    m_params_file << endl <<" ********compute_geno_given_clst_v*********"  << endl;
	    indcnt=0;
	    
	    for(auto  &it:m_geno_given_clst_v)
	    {
		    m_params_file << "IND :" << indcnt << endl;
		    marker =0;
		    
		    for(auto jt : it)
		    {
			  m_params_file << endl<<"Marker "<< marker <<endl;
			  for(auto kt : jt)
			  {
			      m_params_file << " "<<kt << " ";
			  }
			  ++marker;
		    }
		    ++indcnt;
	      }
		m_params_file <<endl << " **********Fwd probs**********"  << endl;
		indcnt=0;  
		for(auto &kvp: m_fwd_probs)
		{
		  m_params_file << "Ind " << indcnt << endl ;
		  marker = 0;
		    for(auto &kv:kvp) 
		    {
			m_params_file << "Marker " << marker <<endl;
			for(auto p:kv)
			  m_params_file << " "<<p <<  " ";  
			++marker;
			m_params_file << endl;
		    }
		  ++indcnt;
		}
	  
		  m_params_file <<endl << " **********Bckwd probs**********"  << endl;
		indcnt=0;  
		for(auto &kvp: m_bckwd_probs)
		{
		  m_params_file << "Ind " << indcnt << endl ;
		  marker = 0;
		    for(auto &kv:kvp) 
		    {
			m_params_file << "Marker " << marker <<endl;
			for(auto p:kv)
			m_params_file << " "<<p <<  " ";  
			++marker;
			m_params_file << endl;
		    }
		  ++indcnt;
		}
	  
		    m_params_file <<endl << " **********compute_clst_given_geno_v**********"  << endl;
		indcnt=0;  
		for(auto &kvp: m_clst_given_geno_v)
		{
		      m_params_file << endl <<"Ind " << indcnt << endl ;
		      ++indcnt;
		      marker = 0;
			for(auto &kv:kvp) 
			{
			   	m_params_file << endl<<"Marker " << marker <<endl;
				for(auto p:kv)
				    m_params_file << " " << p <<  " ";  	      
				++marker;				
			}
		       m_params_file << endl;
		}
}

void em_hmm_genotype::clear_variables()
{
  
  m_hap_data.m_trans_prob_bet_clst.clear();
  m_jump_prob.clear();
  m_geno_given_clst_v.clear();
  m_fwd_probs.clear();
  m_bckwd_probs.clear();
  m_clst_given_geno_v.clear();
  m_clst_given_v.clear();
  //m_ev_jump_given_geno.clear();  
}

void em_hmm_genotype::update_statespace(const gsl_vector *v)
{
    update_param(v,0);
   vector<vector<vector<int>>> stateSpaceUpdated;  
   
   m_part_filter.filter_run(m_clst_given_geno_v,m_current_states,m_stateweights,stateSpaceUpdated);
   m_current_states.clear();	      
   m_current_states = stateSpaceUpdated;
   log_param();
}

gsl_vector* em_hmm_genotype::optimise_param_helper()
{
 //cout<< "CHK em_hmm_genotype::optimise_param_helper() START" << endl; 
      int num_Params;  
      num_Params =  (n_markers* n_clusters * 2) + (n_markers-1) ;//theta, alpha, r .
      
      int paramCount = -1;
      gsl_vector *x = gsl_vector_alloc (num_Params);
      double bounds[2]={0,1}; 
      double dtemp,dtemfinal; 
      //list the initial vales of the parameters
    
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	      for(int jCount=0; jCount < n_clusters ;jCount++)
	      {
		  ++paramCount; 
		  dtemp = m_hap_data.m_theta[iCount][jCount];
		  dtemfinal  = -log((bounds[1]-dtemp)/(dtemp-bounds[0])); 
		  gsl_vector_set (x, paramCount,dtemfinal);
		  
	      }
      }
    
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	      for(int jCount=0; jCount < n_clusters; jCount++)
	      {
		  ++paramCount; 
		  dtemp = m_hap_data.m_alpha[iCount][jCount];
		  dtemfinal  = -log((bounds[1]-dtemp)/(dtemp-bounds[0])); 
		  gsl_vector_set (x, paramCount,dtemfinal);
		
	      }
      }    
    
      for(int iCount=0;iCount < (n_markers-1); iCount++)
      {
	    ++paramCount; 
	    dtemp =m_hap_data.m_recombinations[iCount];
	    dtemfinal  = -log((bounds[1]-dtemp)/(dtemp-bounds[0])); 
	    gsl_vector_set (x, paramCount,dtemfinal);	 
      }
     // cout<< "CHK em_hmm_genotype::optimise_param_helper() END" << endl;
   
   return x;  
} 
 
//p(zim|gim, v*)p(jim|gim, v*)*{log[p(gim|zim, jim, v)]+log[p(zim|v)]+log[p(jim|v)]}
double em_hmm_genotype::func_eval_local(const gsl_vector *v)
{
      // cout<< "CHK em_hmm_genotype::func_eval_local() START" << endl; 
      Double  dClstGgenoV, dJumpGgenoV, dGenoGclustV, dClustGv, dJumpGv=0.0, dMarkValue,dIndValue,functionEval = 0.0;        
      vector<Double> Clst_given_Geno_V_CurMrk, Jmp_given_Geno_V_CurMrk, Geno_given_Clst_V_CurMrk,Clst_given_V;    
      int iscaledTimes = 0;
      int iTemp;
      //double dScaler = pow(10,30);
      double dFinal_functionEval;      
      update_param(v,0);
      //log_param();
      combinable<Double> objVal_mark([]() { return 0.0; });  
      vector<int> cur_markerStates;
      for(int indCnt = 0; indCnt< n_individuals;indCnt++)  
      {
	      dIndValue = 0.0;
	      	      
	      for(int iCurMarker = 0; iCurMarker < n_markers; iCurMarker++)
	      {
		      Clst_given_Geno_V_CurMrk = m_clst_given_geno_v[indCnt][iCurMarker];//p(zim|gim, v*)
		    // Jmp_given_Geno_V_CurMrk  =  m_ev_jump_given_geno[indCnt][iCurMarker];//p(jim|gim, v*)
		      Geno_given_Clst_V_CurMrk = m_geno_given_clst_v[indCnt][iCurMarker];//p(gim|zim, jim, v)
		      Clst_given_V = m_clst_given_v[indCnt][iCurMarker]; // p(zim,jim|v)
		      cur_markerStates = m_current_states[indCnt][iCurMarker];
		     // for(int iStateCount=0; iStateCount < num_states; iStateCount++)
		       parallel_for(int(0), num_states, [this,&objVal_mark,&Geno_given_Clst_V_CurMrk,&cur_markerStates,
				    &Clst_given_V,&Clst_given_Geno_V_CurMrk](int iStateCount)throw()  
		      {
			      //dMarkValue +=  dClstGgenoV*(log(dGenoGclustV)+log(dClustGv));//+log(dJumpGv));
			      objVal_mark.local() +=(Clst_given_Geno_V_CurMrk[iStateCount]*(log(Geno_given_Clst_V_CurMrk[iStateCount]) 
								+ log(Clst_given_V[iStateCount]))/m_total_states[cur_markerStates[iStateCount]].weight);
		      });
		      dIndValue += objVal_mark.combine(plus<Double>());
		      objVal_mark.clear();
		      //dIndValue *= dMarkValue;	    
	      }	  
	      functionEval += dIndValue;
      }     
        iscaledTimes = n_scaler;
        while(iscaledTimes >0)
	{
	    functionEval *= scaling_factor;
	    --iscaledTimes;
	}
	
	functionEval *= -1;
    
      return (double)functionEval;
}

void em_hmm_genotype::diff_eval_local(const gsl_vector *v,gsl_vector *df)
{
 // cout<< "CHK em_hmm_genotype::diff_eval_local() START" << endl; 
  int paramCount = -1;
   //update_param(v,0);
 // log_param();
 // we differentiate function  w.r.t params AT a given marker and given cluster 
      get_theta_gradient(df, &paramCount);//,iCurMarker);
      get_alpha_gradient(df, &paramCount);//,iCurMarker);
      get_recomb_gradient(df,&paramCount);//,iCurMarker);
  //cout<< "CHK em_hmm_genotype::diff_eval_local() END" << endl; 
}

// differentiating wrt to theta_km
// p(zim|gim, v*)p(jim|gim, v*)*{1/p(gim|zim, jim, v)}{d/dtheta_km}[p(gim|zim, jim, v)]
void em_hmm_genotype::get_theta_gradient(gsl_vector *df, int *iParamCount)//, int iCurMarker)
{
  //cout<< "CHK em_hmm_genotype::get_theta_gradient() START" << endl; 
  Double df_theta_km;
 
  vector<int> currentMark_states,tempVec,states_relevant,toStVec, haplo_vector;
  vector<Double> Clst_given_Geno_V_CurMrk, Clst_given_V, Geno_given_Clst_V_CurMrk;
  vector<double> thetaCurMarker;
  vector<vector<int>> permuted_haplotypes;
  Double  dtot_IndValue;
  int toState,loopEnd,iscaledTimes;
  vector< vector<std::pair<int,int>>> Map_Haplos_states; 
    
  combinable<Double> gradtheta_state([]() { return 0.0; });    
  
  for(int iCurMarker = 0; iCurMarker < n_markers; ++iCurMarker)
  {
	    thetaCurMarker = m_hap_data.m_theta[iCurMarker];
	    for(int iCurCluster=0; iCurCluster < n_clusters; ++iCurCluster)
	    {
		      df_theta_km = 0.0;
		      for(int indCnt = 0; indCnt< n_individuals;indCnt++)
		      {
				// constant part of the entire expression		   
				currentMark_states = m_current_states[indCnt][iCurMarker];
				Clst_given_Geno_V_CurMrk = m_clst_given_geno_v[indCnt][iCurMarker];//p(zim|gim, v*)
			    
				Geno_given_Clst_V_CurMrk =  m_geno_given_clst_v[indCnt][iCurMarker];//p(gim|zim, v)
				Clst_given_V = m_clst_given_v[indCnt][indCnt];
				haplo_vector = m_hap_data.m_input.m_genotypes[indCnt][iCurMarker];	   
				dtot_IndValue = 0.0;
			  
				// compute {d/dtheta_km}[p(gim|zim, jim, v)]
				//get all the states with our cluster  'iCurCluster' at the current marker iCurMarker ; rest of them become zero during differentiation.        
				for(int iStateCount=0; iStateCount < num_states; ++iStateCount)
				{
					toState = currentMark_states[iStateCount]; 
					toStVec = m_total_states[toState].values;
					
					if(std::find(toStVec.begin(), toStVec.end(), iCurCluster) != toStVec.end())
					{
					      states_relevant.push_back(iStateCount);
					}
				  }
				  loopEnd = states_relevant.size();
				
				  //we permute the haplotypes,  but not the states
				  vector_permutation(haplo_vector,tempVec,permuted_haplotypes);
			  
				parallel_for(int(0), loopEnd, [this,&Geno_given_Clst_V_CurMrk,&states_relevant,&currentMark_states,&Clst_given_Geno_V_CurMrk,
						&thetaCurMarker, &permuted_haplotypes,&Clst_given_V,&gradtheta_state,&iCurCluster,&indCnt](int stateCount)throw()   
				  {
					  vector<Double> variable_vector_var,variable_vector_const;
					  int icounter,ihaplo;
					  vector< vector<std::pair<int,int>>> Map_Haplos_states;
					  Double const_value,var_value,dtot_stateValue=0.0,dTemphaplo,dThetaOne, dTempvar,dConstantValue=1.0; 
					
					Double dGenoGclustV = Geno_given_Clst_V_CurMrk[states_relevant[stateCount]];					
					vector<int> cur_state = m_total_states[currentMark_states[states_relevant[stateCount]]].values;
					  
					  //Compute the constant value at this individualData
					  for(int markCnt=0; markCnt < n_markers; markCnt++)//p(zim|gim, Jim,v*)
					      dConstantValue *=  Clst_given_Geno_V_CurMrk[states_relevant[stateCount]];
					  //map the state and haplo combi in a unique way,  for eg permuted_haplotypes[i]={0,0,1,1} state={2,3,3,4} map= {(0,2),(0,3),(1,3),(1,4)}
					  get_map_vectors(permuted_haplotypes, cur_state, Map_Haplos_states);
							    
					    for(auto &kvp:Map_Haplos_states)
					    {
							const_value = 1.0;
							var_value = 0.0;
							
							for(auto clst_it:kvp)
							{
								dThetaOne = thetaCurMarker[clst_it.second];
								ihaplo =   clst_it.first;
								// Equation (1) in Scheet & Stephans 2006  article			
								dTemphaplo = pow(dThetaOne,ihaplo);
								dThetaOne = 1.0-dThetaOne;
								ihaplo = 1.0-ihaplo;
								dTemphaplo *= pow(dThetaOne,ihaplo);			  
							    
								if(clst_it.second != iCurCluster)
								{
								      const_value = const_value * dTemphaplo;
								}
								else
								{
								      variable_vector_const.emplace_back(dTemphaplo);
								      if (ihaplo == 0)
									  variable_vector_var.emplace_back(-1);
								      else
									  variable_vector_var.emplace_back(1);
								}
							}
							icounter = 0;		     
							for(auto kp:variable_vector_var)
							{
								dTempvar = kp;
								for (int jcounter = 0; jcounter < variable_vector_const.size();++jcounter)
								{
								      if(jcounter != icounter )
								      {
									  dTempvar = dTempvar * variable_vector_const[jcounter];
								      }      
								}  
								var_value = var_value + dTempvar; 
								++icounter;
							  }
							    
							  variable_vector_var.clear();
							  variable_vector_const.clear();						    
							  //derivative at current permutation of haplo state map 
							  const_value *= Clst_given_V[stateCount];
							  dTempvar = const_value*var_value;	       
							  //sum of all possible permutations 
							  dtot_stateValue +=  dTempvar;	
					    }
					    
					  // include the constant part of the entire expression
					    dtot_stateValue  *=  (dConstantValue/dGenoGclustV); 	
					    // sum of all derivative for all the relevant states
					    gradtheta_state.local() += dtot_stateValue;
				      });
			      
				      df_theta_km +=gradtheta_state.combine(plus<Double>());
				      gradtheta_state.clear();	
				      states_relevant.clear();
				      permuted_haplotypes.clear();
			}//current ind ends here   
			iscaledTimes = n_scaler;
			 while(iscaledTimes>0)// && abs(functionEval)<dThreshold)
			  {
				  df_theta_km *= scaling_factor;
				  --iscaledTimes;	    
			  } 
			gsl_vector_set (df, ++(*iParamCount), (double)df_theta_km);
	    }
  }
      // cout <<  endl;
  // cout<< "CHK em_hmm_genotype::get_theta_gradient() END" << endl; 
}

//differentiating wrt to alpha_km
// p(zim|gim,v*)p(jim|gim,v*)*{1/p(zim|v)}{d/dalpha_km}[p(zim|v))]
void em_hmm_genotype::get_alpha_gradient(gsl_vector *df,int *iParamCount)//,int iCurMarker)
{
   //cout<< "CHK em_hmm_genotype::get_alpha_gradient() START" << endl; 
      
      Double df_alpha_km=0.0,dtot_IndValue;
      int loopEnd,iscaledTimes=0;
      concurrent_vector<int> states_relevant;
      //vector<int> toStVec,states_relevant;
      vector<Double> Clst_given_Geno_V_CurMrk,Clst_given_V;
      combinable<Double> gradalpha_state([]() { return 0.0; });   
      vector<int> cur_markerStates,prevMark_states;
       // for  marker=0
      for(int iCurCluster=0; iCurCluster < n_clusters ;++iCurCluster)
      {
		  df_alpha_km=0.0;
		  for(int indCnt = 0; indCnt< n_individuals;++indCnt)  
		  {	
			      cur_markerStates = m_current_states[indCnt][0];
			      //get all the states with our cluster  'iCurCluster' at the current marker iCurMarker ; rest of them become zero during differentiation.        
			      parallel_for(int(0), num_states, [this,&states_relevant,&iCurCluster,&cur_markerStates](int iStateCount)throw()       
			      {
				      vector<int> toStVec = m_total_states[cur_markerStates[ iStateCount]].values;
				      if(std::find(toStVec.begin(), toStVec.end(), iCurCluster) != toStVec.end())
				      {
					  states_relevant.emplace_back(iStateCount);
				      }
				});
				loopEnd = states_relevant.size();
				Clst_given_V =m_clst_given_v[indCnt][0];
				
				parallel_for(int(0), loopEnd, [this,&iCurCluster,&states_relevant,&indCnt,&gradalpha_state,&Clst_given_V,&cur_markerStates](int stateCount)throw()
				{
					vector<int> cur_stateVec;
					 cur_stateVec = m_total_states[cur_markerStates[states_relevant[stateCount]]].values;
					 Double const_value =1.0, var_value =0.0,dTemp,dConstantValue=1.0;
				
					//Compute the constant value at this individualData
					for(int markCnt=0; markCnt < n_markers; markCnt++)
					{
						dConstantValue *=  m_clst_given_geno_v[indCnt][markCnt][states_relevant[stateCount]];
					}
					 for(auto  &iter:cur_stateVec)
					 {
						if(iter != iCurCluster )
						{
						      const_value = const_value*m_hap_data.m_alpha[0][iter];
						}
						else
						{
							++var_value;
						  }
				          }
					   var_value = var_value*pow(m_hap_data.m_alpha[0][iCurCluster],(var_value-1)); 
					    dTemp = const_value*var_value;// permuted_states.size()
					    //multiply by entire constant value
					    dTemp *= dConstantValue/Clst_given_V[states_relevant[stateCount]];
					    gradalpha_state.local() += dTemp;
				 });
				df_alpha_km  += gradalpha_state.combine(plus<Double>());
				gradalpha_state.clear();
				states_relevant.clear();
		  }
		 iscaledTimes = n_scaler;
		while(iscaledTimes >0)
		{
		      df_alpha_km *= scaling_factor;
		      --iscaledTimes;
	         }
		 gsl_vector_set (df, ++(*iParamCount), (double)df_alpha_km);
		 
         } // marker zero ends here
	
	//differentiate p(zim|v) at m >0,  which is function for Transition probs as Eq 9 in Scheet & Stephans 2006
	for(int iCurMarker = 1; iCurMarker < n_markers; ++iCurMarker)    
	{
		for(int iCurCluster=0; iCurCluster < n_clusters; ++iCurCluster)      
		{
			df_alpha_km =0.0;		    
			for(int indCnt = 0; indCnt< n_individuals;++indCnt)  
			{
			       cur_markerStates = m_current_states[indCnt][iCurMarker];
			       prevMark_states = m_current_states[indCnt][iCurMarker-1];
				//get all the states with our cluster  'iCurCluster' at the current marker iCurMarker ; rest of them become zero during differentiation.        
				parallel_for(int(0), num_states, [this,&states_relevant,&iCurCluster,&cur_markerStates](int iStateCount)throw()       
				{
					vector<int> toStVec = m_total_states[cur_markerStates[ iStateCount]].values;
					if(std::find(toStVec.begin(), toStVec.end(), iCurCluster) != toStVec.end())
					{
					    states_relevant.emplace_back(iStateCount);
					}
				  });
				  loopEnd = states_relevant.size();
				  Clst_given_V =m_clst_given_v[indCnt][iCurMarker];
				  
				 parallel_for(int(0), loopEnd, [this,&iCurCluster,&states_relevant,&indCnt,&gradalpha_state,&Clst_given_V,
					      &cur_markerStates,&prevMark_states,&iCurMarker](int stateCount)throw()
				{
					 vector<int>  cur_stateVec = m_total_states[cur_markerStates[states_relevant[stateCount]]].values;
					 combinable<Double> cur_state([]() { return 0.0; });   
				         Double dcurStateValue =0.0,dTempvar,dconst_sum=0.0,dConstantValue=1.0;
				    
				      
					 parallel_for(int(0), num_states, [this,&cur_stateVec,&prevMark_states,&iCurMarker,&cur_state,&iCurCluster](int fromCount)throw()
					 {
						    vector<int>  tempVec;
						    vector<int> frm_stateVec = m_total_states[prevMark_states[fromCount]].values;
						    vector<vector<int>> permuted_states;
						    vector<Double> variable_vector_const,variable_vector_var;
						    vector< vector<std::pair<int,int>>> Map_fromTo_states; 
						    vector_permutation(frm_stateVec,tempVec,permuted_states);   						    
						    get_map_vectors(permuted_states, cur_stateVec, Map_fromTo_states);
						    Double const_value,var_value,dTemp;
						    int icounter;
						    for(auto &kvp:Map_fromTo_states)
						    {
							    const_value = 1.0;
							    var_value =0.0;
							    
							    for(auto from_to:kvp)
							    {
								  if(from_to.second != iCurCluster)
								  {
									const_value *=  m_hap_data.m_trans_prob_bet_clst[iCurMarker-1][from_to];
								  }
								  else
								  {    
									dTemp = exp(-1*m_hap_data.m_recombinations[iCurMarker-1]*m_hap_data.m_input.m_physical_distances[iCurMarker-1]);
									//dTemp = exp(-1*0.0001*m_hap_data.m_input.m_physical_distances[iCurMarker-1]);
									variable_vector_const.emplace_back(m_hap_data.m_trans_prob_bet_clst[iCurMarker-1][from_to]);
									variable_vector_var.emplace_back(1.0-dTemp);			      
								  }
							    }
							    icounter = 0; 
							  for (auto kv:variable_vector_var)
							  {
								dTemp = kv;
								for (int jcounter = 0;jcounter<variable_vector_const.size(); ++jcounter)
							      { 
								    if (icounter !=jcounter)
									dTemp = dTemp*variable_vector_const[jcounter];
							      }
							      ++icounter;
							      var_value = var_value + dTemp;
							}
							
							variable_vector_var.clear();
							variable_vector_const.clear();
							dTemp = const_value*var_value;
							cur_state.local() += dTemp;	
						    }
					});
					 dcurStateValue = cur_state.combine(plus<Double>());
					 cur_state.clear();
					 
					   // include constatnt term of the entire expression 
					  for(int markCnt=0; markCnt < n_markers; markCnt++)
					  {
					      dConstantValue *=  m_clst_given_geno_v[indCnt][markCnt][states_relevant[stateCount]];	
					  }
					  dConstantValue /=Clst_given_V[states_relevant[stateCount]];
					   dcurStateValue *= dConstantValue;						 
					  gradalpha_state.local() += dcurStateValue;
				});
				 df_alpha_km +=gradalpha_state.combine(plus<Double>());
				 gradalpha_state.clear();	
				 
			}//ind ends
			  iscaledTimes = n_scaler;
			  while(iscaledTimes >0)
			  {
				  df_alpha_km *= scaling_factor;
				--iscaledTimes;
			  }
				
			  gsl_vector_set (df, ++(*iParamCount), (double)df_alpha_km);
			  states_relevant.clear();
			 
		}//cluster ends
	}//marker ends
  
 //cout<< "CHK em_hmm_genotype::get_alpha_gradient() END" << endl; 
}

//differentiating wrt to r_m,
// p(zim|gim, v*)p(jim|gim, v*)*{1/p(zim|v)}{d/dr_m}[p(zim|v))]
void em_hmm_genotype::get_recomb_gradient(gsl_vector *df,int *iParamCount)//,int iCurMarker)
{
  //cout<< "CHK em_hmm_genotype::get_recomb_gradient() START" << endl; 
  
	Double df_r_m=0.0;
	double df_r_m_final;
	vector<double> alphaCurMarker;
	vector<Double> Clst_given_V,Clst_given_Geno_V_CurMrk;
	combinable<Double> gradrecomb_state([]() { return 0.0; });         
	int iscaledTimes;
	vector<int> currentMark_states,prevMark_states;
	for(int iCurMarker = 1;iCurMarker<n_markers;++iCurMarker) 
	{
		alphaCurMarker = m_hap_data.m_alpha[iCurMarker];        
		df_r_m =0.0;   
		for(int indCnt = 0; indCnt< n_individuals;++indCnt)  
		{	
			currentMark_states = m_current_states[indCnt][iCurMarker];
			prevMark_states = m_current_states[indCnt][iCurMarker-1];
			Clst_given_V =m_clst_given_v[indCnt][iCurMarker];
			
			parallel_for(int(0), num_states, [this,&iCurMarker,&alphaCurMarker,&currentMark_states,&prevMark_states,
			&gradrecomb_state,&Clst_given_V,&indCnt](int iStateCount)throw()
			{
				  vector<int> cur_stateVec = m_total_states[currentMark_states[iStateCount]].values;
				  Double dcurStateValue =0.0,dTempk,const_value,dConstantValue=1.0;		     
				  combinable<Double> cur_state([]() { return 0.0; });		
				  
				  parallel_for(int(0), num_states, [this,&iCurMarker,&alphaCurMarker,&prevMark_states,
				  &cur_state,&cur_stateVec](int fromCount)throw()
				  {
					    Double var_value,dTemp;
					    vector<int> frm_stateVec,tempVec;
					    vector<vector<int>> permuted_states;
					    vector< vector<std::pair<int,int>>> Map_fromTo_states;
					    vector<Double> variable_vector_const,variable_vector_var;
					    frm_stateVec = m_total_states[prevMark_states[fromCount]].values;
					    vector_permutation(frm_stateVec,tempVec,permuted_states);   
					    get_map_vectors(permuted_states, cur_stateVec, Map_fromTo_states);
					    int icounter;
						  //differentiating wrt to r_m,  the transition probs as Eq 9 in Scheet & Stephans 2006
					    for(auto &kvp:Map_fromTo_states)
					    {
						    var_value =0.0;
						    for(auto from_to:kvp)
						    {
							  dTemp = exp(-1*m_hap_data.m_recombinations[iCurMarker-1]*m_hap_data.m_input.m_physical_distances[iCurMarker-1]);
							  if(from_to.first == from_to.second) //if to == from 
							  {
								dTemp = dTemp*(alphaCurMarker[from_to.second]-1);
							  }
							  else
							  {    
								dTemp = dTemp*(alphaCurMarker[from_to.second]);
							  }
							  variable_vector_const.emplace_back(m_hap_data.m_trans_prob_bet_clst[iCurMarker-1][from_to]);
							  variable_vector_var.emplace_back(dTemp);		  
						    }
						    icounter = 0; 
						  
						    for (auto kv:variable_vector_var)
						    {
							  dTemp = kv;
							  for(int jcounter = 0;jcounter<variable_vector_const.size(); ++jcounter)
							  { 
								if (icounter !=jcounter)
								{
								      dTemp = dTemp*variable_vector_const[jcounter];
								}
							  }
							  ++icounter;
							  var_value = var_value + dTemp;
						    }
						    variable_vector_var.clear();
						    variable_vector_const.clear();
						    //dcurStateValue = dcurStateValue+var_value;
						    cur_state.local() += var_value;
					    }					  
				  });
				  dcurStateValue = cur_state.combine(plus<Double>());
				  cur_state.clear();
				  //include constatnt term of the entire expression	
				  for(int markCnt=0; markCnt < n_markers; markCnt++)
					dConstantValue *=  m_clst_given_geno_v[indCnt][markCnt][iStateCount];
				  
				  dcurStateValue = dConstantValue / Clst_given_V[iStateCount];
				  
				   gradrecomb_state.local() += dcurStateValue;
			});    
			df_r_m += gradrecomb_state.combine(plus<Double>());
			gradrecomb_state.clear();
		}//ind ends here
		
		iscaledTimes = n_scaler;	      
		while(iscaledTimes >0)
		{
			df_r_m *= scaling_factor;
			--iscaledTimes;
		}
		gsl_vector_set (df, ++(*iParamCount),(double)df_r_m);
	}// marker ends here
  //cout<< "CHK em_hmm_genotype::get_recomb_gradient() END" << endl;  
}

void em_hmm_genotype::update_param(const gsl_vector *x,int status)
{
  //cout<< "CHK em_hmm_genotype::update_param() START" << endl; 
 
    int paramCount = -1;
    double dTemp;
    vector<double> dTempVec;
    m_hap_data.m_theta.clear();
    m_hap_data.m_alpha.clear();
    m_hap_data.m_recombinations.clear();
    double bounds[2]={0,1}; 
    //update the model parameters
    for(int iCount=0;iCount < n_markers; iCount++)     
    {    
	for(int jCount=0; jCount < n_clusters ;jCount++)
	{
	      ++paramCount;
	      dTemp = gsl_vector_get (x, paramCount);	
	      //transform variable to put it inside interval bounds
	      dTemp=bounds[0]+(bounds[1]-bounds[0])/(1+exp(-dTemp));
	      dTempVec.emplace_back( dTemp);		  
	}
	m_hap_data.m_theta.emplace_back(dTempVec);
	dTempVec.clear();
    }
    
    for(int iCount=0;iCount < n_markers; iCount++)     
    {
	  for(int jCount=0; jCount < n_clusters ;jCount++)
	  {
		++paramCount;
		dTemp = gsl_vector_get (x, paramCount);
		//transform variable to put it inside interval bounds
	        dTemp=bounds[0]+(bounds[1]-bounds[0])/(1+exp(-dTemp));
		dTempVec.emplace_back(dTemp);	
	  }
	  m_hap_data.m_alpha.emplace_back(dTempVec);
	  dTempVec.clear();
    }
    
    for(int iCount=0;iCount < (n_markers-1); iCount++)     
    {
	 ++paramCount;
	 dTemp = gsl_vector_get (x, paramCount);
	 //transform variable to put it inside interval bounds
	  dTemp=bounds[0]+(bounds[1]-bounds[0])/(1+exp(-dTemp));
	  m_hap_data.m_recombinations.emplace_back(dTemp);
      } 
  
     // update_statespace();
      update_HMM_param();   
      m_hap_data.log_param_hap();  
    //  log_param();
  
    //cout<< "CHK em_hmm_genotype::update_param() END" << endl; 
}

void em_hmm_genotype::resolve_phase()
{
  Phase PhaseObj;
  vector<vector<vector<int>>> orderedStates;
  vector<vector<vector<int>>> final_Phase;
  
  PhaseObj.resolvePhase(m_fwd_probs,m_clst_given_geno_v,m_hap_data,m_total_states,m_current_states,orderedStates,final_Phase);
  
  log_results(orderedStates,final_Phase);
}

void em_hmm_genotype::log_results(vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)
{
   ofstream clusters_file,phase_file; 
   int ind=0;
   clusters_file.open("Phase_clusters.txt");
   phase_file.open("Phased_haplo.txt");
   vector<vector<int>> indData;

  for(int indCount=0;indCount<n_individuals;++indCount)
   {
      clusters_file << "Ind " << indCount << endl;
      indData = orderedStates[indCount];
      for(int hapCount=0;hapCount<n_ploidy;++hapCount)
      {
	  for(int markCount=0;markCount<n_markers;++markCount)
	  {
	    clusters_file << indData[markCount][hapCount] << " ";
	  }
	  clusters_file << endl;
      }
   }
   for(int indCount=0;indCount<n_individuals;++indCount)
   {
      phase_file << "Ind " << indCount << endl;
      indData = final_Phase[indCount];
      for(int hapCount=0;hapCount<n_ploidy;++hapCount)
      {
	  for(int markCount=0;markCount<n_markers;++markCount)
	  {
	    phase_file << indData[markCount][hapCount] << " ";
	  }
	  phase_file << endl;
      }
   }
   
   phase_file.close();
   clusters_file.close();
}









