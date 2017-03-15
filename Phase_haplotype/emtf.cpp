#include "emtf.h"

CEmtfgeno::CEmtfgeno(string ipFileName,string paramsFileName)
{
     //m_params_file.open ("Test_Phase_param.txt"); 
     //m_hap_data.initialise_param(ipFileName);
     m_hap_data.initialize(ipFileName);
     
     
     if(!paramsFileName.empty())
	initialise_Simmodel_param(paramsFileName,true);
     else
       initialise_param();
     
     //state space
     gen_combi_states(m_total_states); 
     num_states = m_total_states.size();
     n_transitions = num_states/2;//std::min(num_states/3,30);
     //log_param();   
    // cout << "CEmtfgeno::CEmtfgeno"<<endl;
}
 
CEmtfgeno::~CEmtfgeno()
{
    //m_params_file.close();   
}

void CEmtfgeno::initialise_param()
{
    m_hap_data.initialise_param();
}

void CEmtfgeno::populate_transitions()
{
      //cout << "CEmtfgeno::populate_transitions: START " << endl;
      std::random_device rd;
      std::mt19937 gen(rd());
	//std::default_random_engine generator;
	
      std::uniform_int_distribution<int> distribution(0,num_states-1);
      for(int markCnt=0;markCnt<n_markers;++markCnt) 
      {
	    vector<vector<int>> markData;
	    for(int icounter=0; icounter<num_states;++icounter)
	        markData.emplace_back(vector<int>());
	    parallel_for(int (0), num_states,[this,&markData,&distribution,&rd](int fromCount)throw()
	    {
                  int transitionCount=1;
		  int cur_state;
		  //cur_state = distribution(rd);
		  markData[fromCount].emplace_back(fromCount);
		  while(transitionCount<n_transitions)
		  {
		       cur_state = distribution(rd);
		       if(std::find(markData[fromCount].begin(),markData[fromCount].end(),cur_state)
			   ==markData[fromCount].end())
		       {
			   markData[fromCount].emplace_back(cur_state);
			   ++transitionCount;
		       }
		  }	  
	      
	    });
	    m_current_transitions.emplace_back(markData);
      }
   //cout << "CEmtfgeno::populate_transitions: END " << endl; 
}

void CEmtfgeno::initialise_Simmodel_param(string paramsFileName,bool recombDirect)
{
   cout << "CEmtfgeno::initialise_Simmodel_param START" << endl;
   std::ifstream paramsFile(paramsFileName);
   string Strbuffer;
   istringstream buf_stream;
   vector<double> tempData;
   double dTemp;
   //cout << paramsFileName <<endl;
      try
     {
	    if (paramsFile.is_open()) 
	    {
	      
	      getline(paramsFile,Strbuffer);//ignore first line
	       
	       for( int markCount=0;markCount<n_markers;++markCount)
	       {
		  getline(paramsFile,Strbuffer); 		  
		  buf_stream.str(Strbuffer);
		  m_hap_data.m_alpha[markCount] = vector<double>(istream_iterator<double>(buf_stream), istream_iterator<double>()); 
		  buf_stream.clear();
	       }
	       getline(paramsFile,Strbuffer);
	       //getline(paramsFile,Strbuffer);//ignore these lines  
	       for( int markCount=0;markCount<n_markers;++markCount)
	       {
		  getline(paramsFile,Strbuffer);		  
		  buf_stream.str(Strbuffer);
		  m_hap_data.m_theta[markCount] = vector<double>(istream_iterator<double>(buf_stream), istream_iterator<double>()); 
		  buf_stream.clear();
	       }
	       getline(paramsFile,Strbuffer);
	       //getline(paramsFile,Strbuffer);//ignore these lines 
	       getline(paramsFile,Strbuffer);
	       buf_stream.str(Strbuffer);
	       
	       tempData = vector<double>(istream_iterator<double>(buf_stream), istream_iterator<double>()); 
	       if(recombDirect ==true)
	       {
		 for( int markCount=0;markCount<(n_markers-1);++markCount)
		  {
		      m_hap_data.m_recombinations[markCount]=(1-exp(-1* tempData[markCount]*m_hap_data.m_input.m_physical_distances[markCount]));
		  }		   
	      }
	      else
	      {
		  for( int markCount=0;markCount<(n_markers-1);++markCount)
	           {
		      m_hap_data.m_recombinations[markCount]=tempData[markCount];
	           }
	      }
	       
	    }
	    else
	    {
		  cout << "no values in params file, proceeding with random values" << endl;
		  
	    }
     }
    catch(const std::exception& e)
    {
        cout << "Error in reading model_parameters " ;
	std::cout << e.what() << endl;
	exit(1);
    }
   
    cout << "CEmtfgeno::initialise_Simmodel_param END" << endl;
}

void CEmtfgeno::compute_jump_prob()
{
  //cout << "CEmtfgeno::compute_jump_prob() START" <<endl;
        Double jumpProbTemp,jumpProb,dTemp;
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
    //   cout << "CEmtfgeno::compute_jump_prob() END"<< endl;
}
  
void CEmtfgeno::compute_geno_given_clst_v(void)
{
    //cout << "CEmgeno::compute_geno_given_clst_v() START" <<endl;
     vector<int> temp,haplo_vector;
     vector<double> thetaCurMarker;
     vector<vector<int>> genoCurInd, permuted_haplotypes,indStates;
     vector< vector<std::pair<int,int>>> Map_Haplos_states;
    
    m_geno_given_clst_v.clear();
    for(int indCount=0;indCount<n_individuals;++indCount)	
    {
	 m_geno_given_clst_v.emplace_back(vector<vector<Double>>());       
     }
     try
     {
         for(int indCount=0;indCount<n_individuals;++indCount)
	 {
	      genoCurInd = m_hap_data.m_input.m_genotypes[indCount];   
	
	      for(int markCount=0; markCount < n_markers; ++markCount)  
	      { 
		       vector<Double> markerData(num_states);
		       thetaCurMarker = m_hap_data.m_theta[markCount];
		      //get the haplo vector.For eg if geno is 0 for tetraploids,  haplo_vector= {0, 0, 0, 0} as given in the inout file
		      haplo_vector = genoCurInd[markCount];  
		      vector_permutation(haplo_vector,temp,permuted_haplotypes);
		   
		    parallel_for(int(0), num_states, [this,&permuted_haplotypes,&markerData,&thetaCurMarker](int stateCount)throw()
		    {
			  long int ihaplo;
			  Double final_value =0.0,dTemphaplo,perm_values,dThetaOne; 
			  vector<int>  toStatevec;
			  vector< vector<std::pair<int,int>>> Map_Haplos_states;
			  toStatevec = m_total_states[stateCount];
		          get_map_vectors(permuted_haplotypes, toStatevec, Map_Haplos_states);
			  
			  for(auto &it:Map_Haplos_states)
			  {
			      perm_values=1.0;
			       for (int clst_it=0;clst_it<n_ploidy;++clst_it)
			       {
				   dThetaOne = thetaCurMarker[toStatevec[clst_it]];
				   ihaplo =   it[clst_it].first;			
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
  //cout << "CEmgeno::compute_geno_given_clst_v() END" <<endl;
}

//Equation 8  in article
void CEmtfgeno::compute_InitialstateProbs()
{
      m_initial_state_prob.clear();
      m_initial_state_prob = std::vector<Double>(num_states);
      
      parallel_for(int(0), num_states, [this](int stateCount)throw()
     {
	      vector<int> toStVec,tempVec;
	      vector<vector<int>> fromPerms;
	      Double dFinalValue =1.0;
	      toStVec = m_total_states[stateCount];
	      vector_permutation(toStVec,tempVec,fromPerms);
	      for(auto  &iter:toStVec)
	      {
		  dFinalValue = dFinalValue * m_hap_data.m_alpha[0][iter];
	      }
	      dFinalValue = dFinalValue * fromPerms.size();
	      m_initial_state_prob[stateCount]=(dFinalValue);	 				
      });
}
//Equation  9 in article
void CEmtfgeno::compute_clust_given_v()
{
      //cout << "CEmtfgeno::compute_clust_given_v(): START " << endl;
      //m_clst_given_v.clear();
      map< pair<int,int>,Double> transprob_curMarker;  
       try
       { 
		for(int markCnt=0; markCnt < n_markers-1;++markCnt)
		{
		      parallel_for(int(0), num_states, [this,&markCnt](int stateCount)throw()
		      {
				  vector<int> toStVec;			    
				  toStVec  = m_total_states[stateCount];
				  vector<Double> state_data(n_transitions); 
				  Double sum_V, init_prob=0.0;
				 parallel_for(int(0), n_transitions,[this,&toStVec,&state_data,&markCnt,&stateCount](int fromCount)throw()
			        {
				      vector<int> frmStVec = m_total_states[m_current_transitions[markCnt][stateCount][fromCount]];
				      vector<int> tempVec;
				      vector<vector<int>> fromPerms;
				
				       //get the permutations of the from state vector
					vector_permutation(frmStVec,tempVec,fromPerms);
					map<pair<int,int>,Double> markerData;  
					Double transValue = 0.0,tempTrans;
					
					markerData = m_hap_data.m_trans_prob_bet_clst[markCnt];
					vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
					get_map_vectors(fromPerms, toStVec, Map_Valid_stateTuples);
					
					for(auto &kvp:Map_Valid_stateTuples)
					{
						tempTrans = 1.0;
						for(auto kv:kvp)
						{
						    tempTrans *= markerData[kv];
						}
						transValue +=  tempTrans;
					}
					if(transValue < 1e-20)
					  transValue = 1e-20;
					state_data[fromCount] =  transValue;// m_hap_data.compute_trans_prob_bet_clst_tuples(markCnt,fromPerms,toStVec);
					fromPerms.clear();		      
				  });		
				  sum_V = std::accumulate(state_data.begin(),state_data.end(), init_prob);
			         
				  for(auto  &p:state_data)
			             p /= sum_V; 
				  m_clst_given_v[markCnt][stateCount] = state_data;
			 });
	        }	      
      
            }
    catch(const std::exception& e)
    {
	std::cout << e.what() << '\n';
    }
  //cout << "CEmtfgeno::compute_clust_given_v(): END " << endl;
}

 //Equation at the end of Appendix C
void CEmtfgeno::compute_clst_given_geno_v(void)
{
   //cout<< "CHK compute_clst_given_geno_v Start" << endl;  
     vector<Double> FwdProbCurMark, BckwrdProbCurMark;  
      //vector<int> currentMark_states,stateVec,tempVec;
      //vector<vector<int>> perms_frmState,indData;
     // int curState,numPerms;
      Double dTemp,dsumstates;
      
      m_clst_given_geno_v.clear();
      for(int indCount=0;indCount<n_individuals;++indCount)	
     {
	      m_clst_given_geno_v.emplace_back(vector<vector<Double>>());	    
     }
      
      for(int indCount = 0; indCount < n_individuals; ++indCount)
      {
	   for(int markCnt=0; markCnt < n_markers ;++markCnt)
	    {
		    vector<Double> markerData;
		    FwdProbCurMark = m_fwd_probs[indCount][markCnt];
		    BckwrdProbCurMark = m_bckwd_probs[indCount][markCnt];
                    dsumstates =0.0;
		    for(int stateCount=0; stateCount < num_states; ++stateCount)
		    {
			dTemp = FwdProbCurMark[stateCount]*BckwrdProbCurMark[stateCount];
			dsumstates += dTemp;
			markerData.emplace_back(dTemp);
		    }
		    for(int stateCount=0; stateCount < num_states; ++stateCount)
		    {
		       markerData[stateCount] /= dsumstates;
		    }
		    m_clst_given_geno_v[indCount].push_back(markerData);
	      }
     }
}

void CEmtfgeno::compute_fwd_prob(void)
{
  vector<Double> genoGivenClst_cur; 
 // vector<double> AlphaCurMarker;    
  //m_fwd_probs.clear();
 
  for(int indCount=0;indCount<n_individuals;++indCount)
  {
        /*m_fwd_probs.emplace_back(vector<vector<Double>> ());  
	for(int marker=0;marker<n_markers;++marker)
	{
	     m_fwd_probs[indCount].emplace_back(vector<Double>(num_states));
	 }*/
	 genoGivenClst_cur = m_geno_given_clst_v[indCount][0];
	parallel_for(int(0), num_states, [this,&genoGivenClst_cur,&indCount](int stateCount)throw()    //&AlphaCurMarker,
	{
	   vector<int> tostatevec = m_total_states[stateCount];
	   m_fwd_probs[indCount][0][stateCount]=m_initial_state_prob[stateCount]*genoGivenClst_cur[stateCount];
	});
   }
   //vector<Double>FwdProbPrevMarker;
   //for(int indCount=0;indCount<n_individuals;++indCount) 
   for(int markCnt=1; markCnt <n_markers  ;markCnt++)//&AlphaCurMarker,
   {
	//vector<Double>FwdProbPrevMarker,genoGivenClst_cur;	
	parallel_for(int(0), n_individuals, [this,&markCnt](int indCount)throw() 
	{
	     vector<Double> FwdProbPrevMarker = m_fwd_probs[indCount][markCnt-1];
	     vector<Double> genoGivenClst_cur = m_geno_given_clst_v[indCount][markCnt];
	      parallel_for(int(0), num_states, [this,&markCnt,&indCount,&genoGivenClst_cur,&FwdProbPrevMarker](int toStCnt)throw()    //&jumpProbPrevMarker,&AlphaCurMarker,
	      {
		  Double dtempsum=0.0;
		  int state;
		  vector<Double> dclustGv =m_clst_given_v[markCnt-1][toStCnt];
		  vector<int> trans_prevState=m_current_transitions[markCnt-1] [toStCnt];
		  for(int state_id=0;state_id<n_transitions;++state_id)
		  {
		      dtempsum += FwdProbPrevMarker[trans_prevState[state_id]]*dclustGv[state_id];
		  }
		      m_fwd_probs[indCount][markCnt][toStCnt]= genoGivenClst_cur[toStCnt]*(dtempsum);
	      });
	  });
  }
  
 // cout << endl<< "dip_em_hmm_genotype::compute_fwd_prob END" <<endl;
}

void CEmtfgeno::compute_bckwd_prob(void)
{
      //vector<Double>BckwrdProbNxtMarker,geno_given_clst_v;
     //m_bckwd_probs.clear();
     
     /*for(int indCount=0;indCount<n_individuals;++indCount)
	{
	      m_bckwd_probs.emplace_back(vector<vector<Double>> ());      
	      for(int marker=0;marker<n_markers;++marker)
	      {
		  m_bckwd_probs[indCount].emplace_back(vector<Double>(num_states,1.0));
	      }
	}*/
      //vector<Double>BckwrdProbNxtMarker,geno_given_clst_v;
      
        
	for(int markCnt = n_markers-2 ; markCnt >= 0 ;--markCnt)           
       {
	    //vector<Double>BckwrdProbNxtMarker,geno_given_clst_v;
            parallel_for(int(0), n_individuals, [this,&markCnt](int indCount)throw()	   
	   {
	         vector<Double> BckwrdProbNxtMarker =m_bckwd_probs[indCount][markCnt+1];
	          vector<Double> geno_given_clst_v = m_geno_given_clst_v[indCount][markCnt+1];
	      
	         parallel_for(int(0), num_states, [this,&BckwrdProbNxtMarker,
		  &geno_given_clst_v,&markCnt,&indCount](int fromStCnt)throw()    
	         {
		      vector<Double> clust_given_v = m_clst_given_v[markCnt][fromStCnt];
		      Double dtempsum=0.0;  
		   
		      vector<int>trans_nextstate = m_current_transitions[markCnt][fromStCnt];
		      for(int state_id=0;state_id<n_transitions;++state_id)
		      {
			  dtempsum += clust_given_v[state_id]*BckwrdProbNxtMarker[trans_nextstate[state_id]]
			                          *geno_given_clst_v[trans_nextstate[state_id]];
		      }
		      m_bckwd_probs[indCount][markCnt][fromStCnt]= dtempsum;    
		 });
	    });
      }
 }

void CEmtfgeno::log_param()
{
    int marker =0, indcnt=0,state;  
   // ofstream logFile;
    /*if(!fName.empty())
    {
      m_params_file.close();
      m_params_file.open(fName);
    }*/
/*  m_params_file << endl<<"***** current transitions*****"  << endl; 
   for(auto it:m_current_transitions)
   {
       for(auto jt:it)
       {
	 for(auto kt:jt)
	   m_params_file << kt<< " " ;
	 m_params_file << endl;
       }
  }*/
    m_params_file << endl<<"***** Initialstateprobs*****"  << endl;
    for(auto kt : m_initial_state_prob)
	m_params_file << " "  <<kt<< " " ;
   /* m_params_file << endl<<"***** compute_clust_given_v*****";
  
    for(auto  &it:m_clst_given_v)
    { 
	  m_params_file<<endl << "marker: " << marker ;
	  state=0;
	  for(auto jt : it)
	  {   
	   // m_params_file <<endl << "state: " << state << endl ;
	    m_params_file<<endl;
	    for(auto kt : jt)
		  m_params_file << " "  <<kt<< " " ;
	    ++state;
	  }
	  ++marker;
    } */
   m_params_file << endl <<" ********Jump_prob*********";
      
    for(int mCount=0; mCount < n_markers-1; mCount++)
    {
         m_params_file << endl;//<< " marker = "  << mCount<< endl ;

	 for(int jCount=0; jCount <= n_ploidy; jCount++)
	  {
	       m_params_file << m_jump_prob[mCount][jCount] << " ";
	  }
    }
    m_params_file << endl <<" ********compute_geno_given_clst_v*********";
    indcnt=0;
    for(auto  it:m_geno_given_clst_v)
    {
	    m_params_file <<  endl << "IND :" << indcnt;
	    marker =0;
	    for(auto jt : it)
	    {
		  m_params_file << endl;//<< "Marker "<< marker <<endl;
		  for(auto kt : jt)
		  {
		      m_params_file << kt << " ";
		  }
		  ++marker;
	     }
	    ++indcnt;
	   
     }
     m_params_file << endl<<" **********Fwd probs**********";
      indcnt=0;  
      for(auto &kvp: m_fwd_probs)
      {
	    m_params_file << endl << "Ind " << indcnt;
	    marker = 0;
	    for(auto &kv:kvp) 
	    {
		  //m_params_file << endl<< "Marker " << marker <<endl;
	      m_params_file<<endl;
		  for(auto p:kv)
		      m_params_file << p <<  " ";  
		  ++marker;
		 
	      }
	    ++indcnt;
	 }
       
      m_params_file <<endl<<" **********Bckwd probs**********" ;
      indcnt=0;  
      for(auto &kvp: m_bckwd_probs)
      {
	    m_params_file  << endl << "Ind " << indcnt ;
	    marker = 0;
	    
	    for(auto &kv:kvp) 
	    {
		    m_params_file << endl ;//<< "Marker " << marker <<endl;
		    for(auto p:kv)
			m_params_file <<p <<  " ";  
		    ++marker;
		   
	      }
	    ++indcnt;
	 }
    
       m_params_file <<endl << " **********compute_clst_given_geno_v**********"  << endl;
      indcnt=0;  
      for(auto &kvp: m_clst_given_geno_v)
      {
	 m_params_file << "Ind " << indcnt << endl ;
	 marker = 0;
	  for(auto &kv:kvp) 
	  {
	      //m_params_file << "Marker " << marker <<endl;
	      for(auto p:kv)
	       m_params_file << " " << p <<  " ";  	      
	      ++marker;
	      m_params_file << endl;
	  }
	 ++indcnt;
       }
     /*m_params_file <<endl << " **********compute_exp_jump_given_geno_v**********"  << endl;
     indcnt=0;  
     for(auto &kvp: m_ev_jump_given_geno)
      {
	 m_params_file << "Ind " << indcnt << endl ;
	 marker = 0;
	  for(auto &kv:kvp) 
	  {
	      for(auto p:kv)
	       m_params_file << " " << p <<  " ";  	      
	      ++marker;
	      m_params_file << endl;
	  }
	 ++indcnt;
       }*/
  //m_params_file.close(); 
}

void CEmtfgeno::log_model_param(string fName,bool recombDirect)
{
  m_hap_data.log_param_hap(fName,recombDirect);
}
    
void  CEmtfgeno::resolve_phase(bool logging)
{
    cout << "CEmtfgeno::resolve_phase() START" << endl;
    //cout << "m_current_transitions[0][0].size() " << m_current_transitions[0][0].size() << endl;
   // m_hap_data.compute_trans_prob_bet_clst();
   // compute_clust_given_v();
    update_HMM_param();
    //log_param();
    CPhasetf PhaseObj;
    vector<vector<vector<int>>> orderedStates;
    vector<vector<vector<int>>> final_Phase;    
    PhaseObj.resolvePhase(m_fwd_probs,m_clst_given_v,m_current_transitions,m_hap_data,m_total_states,orderedStates,final_Phase,logging);      
    log_results(orderedStates,final_Phase);
      
  cout << "CEmtfgeno::resolve_phase() END" << endl;
}

void CEmtfgeno::log_results(vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)
{
  // cout << "em_hmm_genotype::log_results() START" << endl;
   ofstream clusters_file,phase_file; 
   int ind=0;
   clusters_file.open("Phase_clusters.txt");
   phase_file.open("Phased_haplo.txt");
   vector<vector<int>> indData;
   
   clusters_file << n_individuals << endl;
   phase_file << n_individuals << endl;
   clusters_file << n_markers << endl;
   phase_file << n_markers << endl;
 
   for(int indCount=0;indCount<n_individuals;++indCount)
   {
	clusters_file << "# " << m_hap_data.m_input.m_ind_ID[indCount] << endl;
      
	indData = orderedStates[indCount];
      
	for(int hapCount=0;hapCount<n_ploidy;++hapCount)
	{
	      for(int markCount=0;markCount<n_markers;++markCount)
	      {
		clusters_file << indData[markCount][hapCount]+1 << " ";
	      }
	      clusters_file << endl;
	}
    }
   
    for(int indCount=0;indCount<n_individuals;++indCount)
   {
	  phase_file << "# " << m_hap_data.m_input.m_ind_ID[indCount] << endl;
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
   //std::remove("Phase_param.in");
 //  cout << "em_hmm_genotype::log_results() END" << endl;
}

void CEmtfgeno::update_HMM_param()
{
      try
      {
	  
	    std::vector<std::thread> threads;
	     
	    threads.push_back (std::thread(&CEmtfgeno::compute_jump_prob,this));
	
	    // compute_jump_prob(); 
	     m_hap_data.compute_trans_prob_bet_clst();
	     compute_InitialstateProbs(); 
	    
	    if(threads[0].joinable())threads[0].join();	     
	     
	     threads.emplace_back(std::thread(&CEmtfgeno::compute_clust_given_v,this));
	     //compute_clust_given_v();
	     compute_geno_given_clst_v();
	     for(auto &t:threads)if(t.joinable())t.join();             
	     threads.push_back(std::thread(&CEmtfgeno::compute_bckwd_prob,this));
	    compute_fwd_prob();  	
	    //compute_bckwd_prob(); 
	    for(auto &t:threads)if(t.joinable())t.join();
	     threads.push_back(std::thread(&CEmtfgeno::compute_clst_given_geno_v,this));  
	     //compute_clst_given_geno_v();
	     compute_exp_jump_given_geno_v();
	     for(auto &t:threads) if(t.joinable())t.join();
	         
	}
       catch(const std::exception& e)
      {
	    cout << "EXCEPTION CEmtfgeno::update_HMM_param() " << endl; 
	    std::cout << e.what() << '\n';
      }
	    
}

double CEmtfgeno::run_em(int iterations)
{
      cout << "CEmtfgeno::run_em() START" << endl;
      int i = 0;
      bool bIsconverged = false,bTheta=false,bAlpha=false,bRecomb=false;  
      vector<vector<double>> alphaUpdated;
      vector<vector<double>> thetaUpdated;
      vector<double> recombinationsUpdated;
      vector<vector<int>> stateSpaceUpdated;
      //vector<double> alpha_firstMark = m_hap_data.m_alpha[0];  
	    
      string fName = "model_param.txt";
      ofstream myfile(fName);  
      myfile<< " EM algo start:"<<endl;
      
      double dcurrent_expectation;
      double dprev_expectation= std::numeric_limits<double>::min();
      
      try
     {
	    populate_transitions();
	    
	    for(int indCount=0;indCount<n_individuals;++indCount)
	    {
		  m_fwd_probs.emplace_back(vector<vector<Double>> ());  
		  m_bckwd_probs.emplace_back(vector<vector<Double>> ());      
		  m_ev_jump_given_geno.emplace_back(vector<vector<Double>> ()); 
		  for(int marker=0;marker<n_markers-1;++marker)
		  {
		      m_fwd_probs[indCount].emplace_back(vector<Double>(num_states));
		      m_bckwd_probs[indCount].emplace_back(vector<Double>(num_states,1.0));
		      m_ev_jump_given_geno[indCount].emplace_back(vector<Double>(n_clusters));
		  }
		  m_fwd_probs[indCount].emplace_back(vector<Double>(num_states));
		  m_bckwd_probs[indCount].emplace_back(vector<Double>(num_states,1.0));
	    }
	      
	    for(int marker=0;marker<n_markers-1;++marker)
	    {
		vector<vector<Double>> temp;
		for(int state=0; state<num_states;++state)
		    temp.emplace_back(vector<Double>());
		m_clst_given_v .emplace_back(temp);
	    }    
             alphaUpdated.emplace_back(m_hap_data.m_alpha[0]);
             for(int mark=0;mark<n_markers-1;++mark)
	     {
		     alphaUpdated.emplace_back(vector<double>(n_clusters));
	      }
	
	      recombinationsUpdated = vector<double>(n_markers-1);
	      for(int mark=0;mark<n_markers;++mark)
		 thetaUpdated.emplace_back(vector<double>(n_clusters)); 
	   do
	   {
		//compute expecation
	          ++i;
		  update_HMM_param();
		  //log_param();	   
	          
		dcurrent_expectation = compute_expectation();     	    
		  if(std::isnan(dcurrent_expectation))
                  {
			 cout << "dcurrent_expectation is nan" << endl;
			 log_model_param(fName,false);
			 break;
		  }
		
		  bIsconverged = test_likelihood(dprev_expectation,dcurrent_expectation);
		  dprev_expectation=dcurrent_expectation;
		  
		  myfile << "Iter "<< i << " " << dcurrent_expectation <<" converged "<<bIsconverged << endl;      
		  //myfile.close();
		  
		  //Maximize
		  std::vector<std::thread> threads;
		  threads.push_back(std::thread(&CEmtfgeno::filter_transitions,this,i));		  
		  threads.push_back(std::thread(&CEmtfgeno::compute_theta_updated,this,std::ref(thetaUpdated)));
		
                  //compute_theta_updated(thetaUpdated);
		  compute_alpha_updated(alphaUpdated);
		  compute_recombination_updated(recombinationsUpdated);    
		  
		  for(auto &t:threads)
		  if(t.joinable())t.join();       
		  //update 		  
		  m_hap_data.m_theta = thetaUpdated;
		  m_hap_data.m_alpha = alphaUpdated;
		  m_hap_data.m_recombinations = recombinationsUpdated;
		  
		  //log_model_param(fName);
		  //myfile.open(fName,ios::out | ios::app);  	  
		  
	}while(!bIsconverged && i<iterations);
 
  }
  catch(const std::exception&e)
  {
    std::cout << e.what() << '\n';
  }
  myfile.close();
 // log_model_param(fName,true);
  cout << "CEmtfgeno::run_em() END" << endl;
  return dcurrent_expectation;
 // cout << "CEmtfgeno::run_em() END" << endl;
}

void CEmtfgeno::compute_exp_jump_given_geno_v(void)
{
  // cout << "CEmtfgeno::compute_exp_jump_given_geno_v() START" << endl;
     
     vector<Double> BckwrdProbCurMark,Geno_given_Clst_V_CurMrk;
     
      combinable<Double> expJumpind_clust([]() { return 0.0; });  
      Double dFwdProbPrevMark,dtotprobInd,dsummark,dInit=0.0;
       
      for(int indCnt=0;indCnt<n_individuals;++indCnt)
      {
	//  vector<vector<Double>> indData;
	  dtotprobInd = std::accumulate(m_fwd_probs[indCnt][(n_markers-1)].begin(),m_fwd_probs[indCnt][(n_markers-1)].end(),dInit);
	 //our actual iCurMarker is iCurMarker+1 as there are no JUMPS at  iCurMarker=0
	  for(int iCurMarker=0;iCurMarker<(n_markers-1);++iCurMarker)
	  {
	        //vector<Double> markData(n_clusters);
	        Geno_given_Clst_V_CurMrk=m_geno_given_clst_v[indCnt][iCurMarker+1];
	        BckwrdProbCurMark = m_bckwd_probs[indCnt][iCurMarker+1];
	        dFwdProbPrevMark = std::accumulate(m_fwd_probs[indCnt][iCurMarker].begin(),m_fwd_probs[indCnt][iCurMarker].end(),dInit);	     
	        
		//dsumExpecationscurMark=0.0;
	         parallel_for(int(0),n_clusters,[this,&dFwdProbPrevMark,&dtotprobInd,&Geno_given_Clst_V_CurMrk,
		  &BckwrdProbCurMark,&iCurMarker,&indCnt,&expJumpind_clust](int iCurCluster)throw()//
	         {
		      vector<int> states_relevant,curstateVec;
		      int loopEnd,curstateId,iClust,iClusterInd,iJumpCheck;
		      Double dalphaValue,djump=0.0,dtemp,dexp_probcurclst;
		      dexp_probcurclst =0.0;
		      for(int njumps=1;njumps<=n_ploidy;++njumps)
		      {
			     for(int iStateCount=0;iStateCount<num_states;++iStateCount)
		             {
		                    //std::count (myvector.begin(), myvector.end(), 20);
		                    if(std::count(m_total_states[iStateCount].begin(), 
			             m_total_states[iStateCount].end(), iCurCluster) >=njumps)
			            {
			                states_relevant.emplace_back(iStateCount);			      
			            }			       
		             }	
		             loopEnd = states_relevant.size();
			     
			     for(int statecount=0;statecount<loopEnd; ++statecount)
			     {
			            curstateId = states_relevant[statecount];
			            curstateVec = m_total_states[curstateId];
				    djump =0.0;
				    // num of possible jumps can be more than njumps, but only njumps at iCurCluster
				    for(int posjumps=njumps;posjumps<n_ploidy;++posjumps)
				    {   
				            //jump logic
				             iClust = njumps;
				            dalphaValue =pow(m_hap_data.m_alpha[iCurMarker+1][iCurCluster],njumps);
					     iClusterInd=0;
					     iJumpCheck = njumps;
					      while(iClust<posjumps)
					      {
						  if(curstateVec[iClusterInd] != iCurCluster)
						  {
						         dalphaValue *= m_hap_data.m_alpha[iCurMarker+1][curstateVec[iClusterInd]];
							 ++iJumpCheck;
						  }
						 
						 ++iClusterInd;
						 ++iClust;
					      }
					      if(iJumpCheck ==posjumps)
						  djump += m_jump_prob[iCurMarker][posjumps] * dalphaValue;
				    }
				    dexp_probcurclst +=  dFwdProbPrevMark*djump*BckwrdProbCurMark[curstateId] 
				                                   *Geno_given_Clst_V_CurMrk[curstateId];
			     }
		              states_relevant.clear();        			      
			      
		    } 
		    dexp_probcurclst = (dexp_probcurclst+0.00001)/(dtotprobInd+0.00002);
		    m_ev_jump_given_geno[indCnt][iCurMarker][iCurCluster]= dexp_probcurclst;	
		    //expJumpind_clust.local() += dexp_probcurclst;			      
	     });
	     
	       //dsummark = expJumpind_clust.combine(plus<Double>());
	    // expJumpind_clust.clear();
	    // for(auto &it:markData)
		//it /= dsummark;
	      //indData.emplace_back(markData);
	      
	  }
	//  m_ev_jump_given_geno.emplace_back(indData);
      }
    
   
  // dead code!	    
 //cout << "CEmtfgeno::compute_exp_jump_given_geno_v() END" << endl;
}

bool CEmtfgeno::test_likelihood(double old_likelihood,double new_likelihood)
{
 //cout << "CEmtfgeno::test_convergence() Start" << endl;
  bool bIsConverged = false;
  double tolerance=0.0000001;
   if(is_dclose(old_likelihood,new_likelihood))// || (old_likelihood-new_likelihood)>5.0)
   {
       bIsConverged =true;    
   }
 
 // cout << "CEmtfgeno::test_convergence() End " << bIsConverged << endl;  
  return bIsConverged;
}

void CEmtfgeno::compute_theta_updated(vector<vector<double> > &T)
{
  //cout << "CEmtfgeno::compute_theta_updated() Start" << endl; 
 
     //vector<double> thetaCurMarker;
  
     //for(int iCurMarker=0;iCurMarker<n_markers;++iCurMarker)
       parallel_for(int(0),n_markers,[this,&T](int iCurMarker)throw()
     { 
            vector<double> thetaCurMarker = m_hap_data.m_theta[iCurMarker];
	    
	    //for(int iCurCluster=0;iCurCluster<n_clusters;++iCurCluster)        	    
	    parallel_for(int(0),n_clusters,[this,&T,&thetaCurMarker,&iCurMarker](int iCurCluster)throw()         
	    {    
	          vector<int> states_relevant,curStateVec,genoVec,tempVec,tempVec_NotcurClust;
		  vector<Double> clstGivenGenoCurMark;
		  vector<vector<int>> perm_States,perm_NotcurClust,states_dTempNumNum;
		  int loopEnd,genosumInd,icounter,NumOccur,curState,hapcount;
		  Double dtotprobAllInd=0.0,dTempNumterm2,dTempNumterm1,dprod,
		  dTempNumNum,dTempNumdenom,dNumerator=0.0,dtheta_km,dTempNumdenomTemp,genosumcount;
		  dtotprobAllInd=0.0;
		  dNumerator=0.0;
		  states_relevant.clear();
		  
		  for(int iStateCount=0;iStateCount<num_states;++iStateCount)
		  {
		         if(std::find(m_total_states[iStateCount].begin(), 
			  m_total_states[iStateCount].end(), iCurCluster) != m_total_states[iStateCount].end())
			 {
			      states_relevant.emplace_back(iStateCount);			      
			  }
		   }
		   loopEnd = states_relevant.size();	   
		      
		   for(int indCnt=0;indCnt<n_individuals;++indCnt)
		   {
			    clstGivenGenoCurMark = m_clst_given_geno_v[indCnt][iCurMarker];
			    genoVec = m_hap_data.m_input.m_genotypes[indCnt][iCurMarker];
			    
			    //get the gentype as the sum of alleles
			    genosumInd = std::accumulate(genoVec.begin(),genoVec.end(),0);
			    for(int iStateCount=0;iStateCount<loopEnd;++iStateCount)
			    {				      
				      curState = states_relevant[iStateCount]; 
				      curStateVec =  m_total_states[curState];  
				      NumOccur = std::count(curStateVec.begin(),curStateVec.end(), iCurCluster);
				  
				      dTempNumterm2 = clstGivenGenoCurMark[curState]*NumOccur ;
				    
				      dtotprobAllInd +=  dTempNumterm2 ;
				      dTempNumterm1=1.0;
				      if(genosumInd >0 && genosumInd <n_ploidy)
				      {
					      perm_States.clear();
					      vector_permutation(curStateVec,tempVec,perm_States);
					      dTempNumNum = 1.0;
					      dTempNumdenom = 0.0;
					      genosumcount=0;
					      icounter=0;
					      tempVec.clear();
					      hapcount = 0;
					      tempVec_NotcurClust.clear();
					     while(hapcount<n_ploidy)
					      {
						    //read the clusters other than current clusters into another temp vec, permute them and then use them to
						    //to get states beginning with current cluster and then proceed with numerator							
						      if(curStateVec[hapcount] != iCurCluster)
							tempVec_NotcurClust.emplace_back(curStateVec[hapcount]);
						    ++hapcount;
					      }
					      states_dTempNumNum.clear();
					      vector_permutation(tempVec_NotcurClust,tempVec,states_dTempNumNum);
					      
					      for(auto &it:states_dTempNumNum)
					      {
						   it.insert(it.begin(),NumOccur,iCurCluster);    
					      }
					      dTempNumNum =0.0;
					      for(auto &it:states_dTempNumNum)
					      {
						     dprod=1.0;
						     hapcount=0;
						     genosumcount=0;
					             while(hapcount<n_ploidy)
						    {
							  if(genosumcount < genosumInd )
							  {
								  dprod*= thetaCurMarker[it[hapcount]];
								  ++genosumcount;							    
							  }					      
							  else 
							  {
								  dprod*= (1-thetaCurMarker[it[hapcount]]);
							  } 
							  ++hapcount; 
						    }
						    dTempNumNum += dprod;
					      }
					      for(auto &stateVariant:perm_States)
					      {
						    dprod = 1.0;
						    hapcount = 0;
						    genosumcount=0;
						    for(auto &stVarCnt:stateVariant)
						    {
							    if(genosumcount < genosumInd)
							    {
								dprod*= thetaCurMarker[stVarCnt];	
								++genosumcount;
							    }					      
							    else 
							    {
								dprod*= (1-thetaCurMarker[stVarCnt]);								
							    } 
							    ++hapcount; 
						      }
						      dTempNumdenom +=  dprod;
					      }
					      dTempNumterm1 = (dTempNumNum/dTempNumdenom);
					}
					else if(genosumInd == 0)
					  dTempNumterm1 =0.0;
					else
					  dTempNumterm1 =1.0;
					
				dNumerator +=  (dTempNumterm1*dTempNumterm2);
				
				}//state loop ends here
		      }////ind loop ends here
		      
		      dtheta_km =(dNumerator)/(dtotprobAllInd);
		      if(dtheta_km <0.0001 )
			dtheta_km =0.0001;
		      else if(dtheta_km>0.999)
			dtheta_km =0.999;
		      T[iCurMarker][iCurCluster] = (double)(dtheta_km) ;
	      });// cluster loop ends here
	  });
	  
  //cout << "CEmtfgeno::compute_theta_updated() End" << endl;
  
}

void CEmtfgeno::compute_alpha_updated(vector<vector<double> > &T)
{
    parallel_for(int(0),n_markers-1,[this,&T] (int iCurMarker) throw()
   { 
	  Double dtotExpValueatkm=0.0, dInit=0.0;
          double dsummark;
          for(int indCnt=0;indCnt<n_individuals;++indCnt)
	  {
		dtotExpValueatkm += std::accumulate(m_ev_jump_given_geno[indCnt][iCurMarker].begin(),m_ev_jump_given_geno[indCnt][iCurMarker].end(),dInit);
	  }
	  
	  parallel_for(int(0),n_clusters,[this,&T,&dtotExpValueatkm,&iCurMarker] (int iCurCluster)throw()
	  {
	      Double dalpha_km =0.0;
	      for(int indCnt=0;indCnt<n_individuals;++indCnt)
	      {
		dalpha_km += m_ev_jump_given_geno[indCnt][iCurMarker][iCurCluster];
	      }
	      
	      dalpha_km = ((dalpha_km+0.00001)/(dtotExpValueatkm+0.00002));
	      
	      T[iCurMarker+1][iCurCluster]=(double)(dalpha_km);;
	   
           });
            dsummark = std::accumulate(T[iCurMarker+1].begin(),T[iCurMarker+1].end(),0.0);
	   /* for(int iCurCluster=0;iCurCluster<n_clusters;++iCurCluster)
	    {
	        T[iCurMarker+1][iCurCluster] /= dsummark;
	    }		*/

  });    
  //cout << "CEmtfgeno::compute_alpha_updated() End" << endl;
}

void CEmtfgeno::compute_recombination_updated(vector<double> &T)
{
  //cout << "CEmtfgeno::compute_recombination_updated() START" << endl;
  
     parallel_for(int(0),n_markers-1,[this,&T] (int iCurMarker)throw()
    {
	  Double d_rm,dtotExpValueatkm=0.0;
          Double dInit=0.0;
          double dFinalValue;
    
	  for(int indCnt=0;indCnt<n_individuals;++indCnt)
	  {
	      dtotExpValueatkm += std::accumulate(m_ev_jump_given_geno[indCnt][iCurMarker].begin(),m_ev_jump_given_geno[indCnt][iCurMarker].end(),dInit);
	  }  
	      
	  d_rm =  (dtotExpValueatkm)/(n_ploidy*n_individuals*n_clusters);
	  dFinalValue = (double)d_rm;
	  T[iCurMarker]=(dFinalValue);  
    });
  //cout << "CEmtfgeno::compute_recombination_updated() End" << endl; 
}

double CEmtfgeno::compute_expectation()
{
  //cout << "CEmtfgeno::compute_expectation() START" << endl;
  Double dexp_valueAllInd=0.0,dInit=0.0,dIndexpval;
  double dexp_valueFinal;
  
  for(int indCnt=0;indCnt<n_individuals;++indCnt)
  {    
     dexp_valueAllInd += log(std::accumulate(m_fwd_probs[indCnt][(n_markers-1)].begin(),m_fwd_probs[indCnt][(n_markers-1)].end(),dInit));  
     
  }
   
   dexp_valueFinal= (double)dexp_valueAllInd;
  // cout << "CEmtfgeno::compute_expectation() END" << endl;
  return dexp_valueFinal;
}

void  CEmtfgeno::util_hmm_fwd(int marker,int state,int njumps,vector<vector<int>> &jumpData, 
			  vector<Double>&alpha_data)
{
 // cout<< "CHK em_hmm_genotype::util_hmm_fwd() START" << endl;
  vector<double> AlphaMarker;
  int state_jump;
  Double dalpha_temp;
  vector<int> input_comb,copy_Curstate;
   
   
   // WE FIRST work to obtain the jumpData
   vector<int> curState= m_total_states[state];
   vector<int> gotTemp,currentMark_trans;  
   vector< vector<int>> T,Tfinal;
   vector< vector<int>> ::iterator itTemp;
   // fetch all the possible combinations of clusters with given number of jumps in current state
   vector_combination(curState,gotTemp, 0, njumps,T);

   for(int i=0;i<n_clusters;++i)
      input_comb.emplace_back(i);
   vector<int>::iterator pos;   
   
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
		    itTemp = std::find( m_total_states.begin(), m_total_states.end(), iter);
		    
		    if(itTemp!=m_total_states.end())
		    {
			state_jump = std::distance( m_total_states.begin(),itTemp);
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

void CEmtfgeno::filter_transitions(int curIteration)
{
   //cout << "CEmtfgeno::filter_transitions() START " << endl;
   vector<vector<Double>> temp_expTrans;
   
   int n_transitionsPrev=n_transitions;
   if(curIteration >n_iterations/7 && curIteration<n_iterations/3)
     n_transitions = num_states/3; //std::min(num_states/4,20);
   
   if(curIteration >n_iterations/3)
     n_transitions = num_states/4;//std::min(num_states/5,12);
  //  cout << curIteration<< " "<<n_transitions << endl;
   
   for(auto &i:m_current_transitions)
     for(auto &j:i)
   j.erase(j.begin()+n_transitions,j.end());

    //cout << "m_current_transitions[0][0].size() " << m_current_transitions[0][0].size()<< endl;
    int states_eraseAdd  = n_transitions/3;
    Double dInitValue=0.0;
    map< pair<int,int>,Double> transprob_curMarker;
    vector<Double> IndValues(n_individuals);
    vector<vector<Double>> IndTransValues;//(num_states);
    
    //total value of p(gi|v)
    for(int indCnt=0;indCnt<n_individuals;++indCnt)
	  IndValues .emplace_back(std::accumulate(m_fwd_probs[indCnt][n_markers-1].begin(),m_fwd_probs[indCnt][n_markers-1].end(),dInitValue)); 
	
   for(int markCnt=0;markCnt<n_markers-1;++markCnt)
   {
       transprob_curMarker = m_hap_data.m_trans_prob_bet_clst[markCnt];
         temp_expTrans.clear();
	 for(int it=0;it<num_states;++it)
	      temp_expTrans.emplace_back(vector<Double>());
	  IndTransValues.clear();
	  for(int it=0;it<n_individuals;++it)
	  {
	    IndTransValues.emplace_back(vector<Double>(num_states));
	    
	  }
	   

	   // compute temp values required to compute expected transition Values of current transitions common for all states
	   parallel_for(int(0),n_individuals,[this,&markCnt,&IndValues,&IndTransValues](int indCnt)throw()
          {
	   parallel_for(int(0),num_states,[this,&indCnt,&markCnt,&IndTransValues,&IndValues](int stateCount)throw()
	   {
	       IndTransValues [indCnt][stateCount] =(( m_fwd_probs[indCnt][markCnt][stateCount]*m_geno_given_clst_v[indCnt][markCnt+1][stateCount]
		                           *m_bckwd_probs[indCnt][markCnt+1][stateCount])/IndValues[indCnt]); 
	  });
	  
	  });
	   
	   //compute expected transitions
          parallel_for(int(0),num_states,[this,&temp_expTrans,&markCnt,&IndValues,&IndTransValues](int stateCount)throw()
	  {
	        vector<Double> tempTrans(n_transitions);
		
		parallel_for(int(0),n_transitions,[this,&tempTrans,&markCnt,&stateCount,&IndValues,&IndTransValues](int transCount)throw()
		{
		    //Double dInitValue=0.0,finValue=0.0;
		    //Double IndValue;
		    Double transValue = m_clst_given_v[markCnt][stateCount][transCount];
		   combinable<Double> finValue([]() { return 0.0; });  
		    //for(int indCnt=0;indCnt<n_individuals;++indCnt)
		    parallel_for(int(0),n_individuals,[this,&finValue,&transValue,&markCnt,&stateCount,&IndValues,&IndTransValues](int indCnt)throw()
		    {
		         
		         /*finValue.local() +=(( m_fwd_probs[indCnt][markCnt][stateCount]*transValue*m_geno_given_clst_v[indCnt][markCnt+1][stateCount]
		                           *m_bckwd_probs[indCnt][markCnt+1][stateCount])/IndValues[indCnt]);*/ 
			 finValue.local() +=IndTransValues[indCnt][stateCount] *transValue;
		    });        
		    
	             tempTrans[transCount] = finValue.combine(plus<Double>());
		     finValue.clear();
		});
		temp_expTrans[stateCount]= tempTrans;
	  });
	  //sort the vector in descending order and store the sorted indices in transition_indices
	  // sample for new transitions in updated_transitions
	parallel_for(int(0),num_states,[this,&temp_expTrans,&markCnt,&states_eraseAdd,&transprob_curMarker,&n_transitionsPrev](int stateCount)throw()
	  {
	       vector<int>updated_transitions;
                std::vector<int> transition_indices(n_transitions);
	        std::size_t n(0);
		int state_count,state,loopend,state_Found;
		vector<int>::iterator it_index;
	        std::generate(std::begin(transition_indices), std::end(transition_indices), [&]{ return n++; });
	        std::sort(  std::begin(transition_indices), 
			  std::end(transition_indices),
			  [&](int i1, int i2) { return temp_expTrans[stateCount][i1] > temp_expTrans[stateCount][i2]; } );		
		state_count =  states_eraseAdd;
		vector<int>::iterator it=transition_indices.begin();
		loopend = n_transitions- states_eraseAdd;
		vector<int> frmStVec,toStVec;
		toStVec = m_total_states[stateCount];
		Double transition_newCost;
		// look for transitions which are close to the best transititions and then compare their 
		//transition values to those transitions which are being replaced
		for(int loop_count=0;loop_count<loopend ;++it,++loop_count)
		{
		      state = m_current_transitions[markCnt][stateCount][*it];
		      state_Found =0;
		      if((state+1) <= (num_states-1))
		      {
		         it_index = std::find(m_current_transitions[markCnt][stateCount].begin(),m_current_transitions[markCnt][stateCount].end(),(state+1));
			  if(it_index == m_current_transitions[markCnt][stateCount].end())		      
			    if(std::find(updated_transitions.begin(),updated_transitions.end(),(state+1))==updated_transitions.end())
			    {
			        state = state +1;
				state_Found =1;
			    }
		      }
		      if((state-1)>0 && state_Found==0)
		       {
			    it_index = std::find(m_current_transitions[markCnt][stateCount].begin(),m_current_transitions[markCnt][stateCount].end(),(state-1));
		            if(it_index == m_current_transitions[markCnt][stateCount].end())
		               if(std::find(updated_transitions.begin(),updated_transitions.end(),state-1)==updated_transitions.end())
			      {
			        state = state -1;
				state_Found =1;
			      }
		        }
		        if(state_Found ==1)
			{
				    frmStVec = m_total_states[state];    
				    // compute the transition cost of this new state
				    if(frmStVec[0] !=frmStVec[1] && toStVec[0]!=toStVec[1])
				   {
				        transition_newCost=(transprob_curMarker[std::make_pair(frmStVec[0],toStVec[0])]*transprob_curMarker[std::make_pair(frmStVec[1],toStVec[1])])
				        +(transprob_curMarker[std::make_pair(frmStVec[0],toStVec[1])]*transprob_curMarker[std::make_pair(frmStVec[1],toStVec[0])]);
				    }
				     else
				     {
				           transition_newCost = (transprob_curMarker[std::make_pair (frmStVec[0],toStVec[0])]
					      *transprob_curMarker[std::make_pair (frmStVec[1],toStVec[1])]);
				      }
				    //compare with transition cost of state being replaced 
				  if(m_clst_given_v[markCnt][stateCount][loop_count+states_eraseAdd] < transition_newCost)
				  {
					updated_transitions.emplace_back(state);
				   }
			    }
                }
		       	    
		//erase bad states 
		
		transition_indices.erase (transition_indices.begin(),(transition_indices.begin()+n_transitions-updated_transitions.size()));
		std::sort(transition_indices.begin(), transition_indices.end(), std::greater<int>());
	      for(int i:transition_indices)
	      {
		//m_current_transitions[markCnt][stateCount].erase(m_current_transitions[markCnt][stateCount].begin()+n_transitions-1,m_current_transitions[markCnt][stateCount].end());
		m_current_transitions[markCnt][stateCount].erase( std::next( m_current_transitions[markCnt][stateCount].begin(), i ) );		
	      }
		m_current_transitions[markCnt][stateCount].insert(m_current_transitions[markCnt][stateCount].end(),updated_transitions.begin(),updated_transitions.end()); 
                
	  });
	
   }
 //  cout << "m_current_transitions[0][0].size() " << m_current_transitions[0][0].size()<< endl;
   ////cout << "CEmtfgeno::filter_transitions() END " << endl;
}

//Test code
int main(int argc, char *argv[])
{
    string ipFilename,parFName;
    string runsFileName= "runs_info.txt";
    if ( argc < 2 ) 
    {
	// argv[0] is the program name
	//cout<<"usage: "<< argv[0] <<" <input data filename>   <Sim Mod paramsfilename >(optional)\n";
	exit(1);
    }
    else 
    {
	    // argv[1] is  input filename to open
	    ifstream the_file ( argv[1] );
	    
	    if ( !the_file.is_open() )
	    {
		  cout<<"Could not open file\n";
		  exit(1);
	    }
	    else 
	    {
		  ipFilename =  argv[1];
		  the_file.close();		  
	    }
	   if(argc == 3)
           	parFName= argv[2];
	   
    }
 
  CEmtfgeno HMMObj(ipFilename,parFName);
  ofstream O_file(runsFileName);
  ifstream I_file;
  int i_runcounter=0;
  double dexpectation;
  vector<double> expectations;
  string paramsFileName="tempInitialparams.txt";
  tbb::tick_count before = tbb::tick_count::now();
  dexpectation = HMMObj.run_em(n_iterations);//n_iterations
  tbb::tick_count after = tbb::tick_count::now();
 cout << "time taken " << (after - before).seconds() << " seconds" << endl;
  
 //before = tbb::tick_count::now();
  HMMObj.resolve_phase(true);
  //after = tbb::tick_count::now();
  //cout << "time taken resolve_phase " << (after - before).seconds() << " seconds" << endl;
  
  
  /*if(parFName.empty())
  {
	while(i_runcounter<n_runs)
	{
	    ++i_runcounter;
	    dexpectation = HMMObj.run_em(3);
	    expectations.emplace_back(dexpectation);
	    O_file << "run "<< i_runcounter<< endl; 
	    O_file.close();
	    HMMObj.log_model_param(runsFileName,false); 
	    HMMObj.initialise_param();
	    O_file.open(runsFileName,ios::out | ios::app); 
	}
	O_file.close();
	
	i_runcounter =  std::distance(expectations.begin(), std::max_element(expectations.begin(), expectations.end())) << '\n';
	
	//read params from runs_info.txt for the run number in i_runcounter in store in parmsFilename
	I_file.open("runs_info.txt");
	I_file.close();

	HMMObj.initialise_Simmodel_param(paramsFileName,false);
	
    }
   else
   {
      dexpectation = HMMObj.run_em(1);
  }
  tbb::tick_count before = tbb::tick_count::now();
  dexpectation = HMMObj.run_em(n_iterations);
 tbb::tick_count after = tbb::tick_count::now();
 cout << "time taken " << (after - before).seconds() << " seconds" << endl;*
  HMMObj.resolve_phase(true);*/
   
}
