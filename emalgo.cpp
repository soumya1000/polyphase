#include "emalgo.h"


//*************************************************************//
// CLASS EM_HMM_GENOTYPE
//*************************************************************//

em_hmm_genotype::em_hmm_genotype()
{
      //cout << "em_hmm_genotype::em_hmm_genotype() START" << endl;
      m_hap_data.restore_prevstate();
      get_total_stateSpace();  
      // cout << "em_hmm_genotype::em_hmm_genotype() END" << endl;
 }

em_hmm_genotype::~em_hmm_genotype()
{
   //cout << "em_hmm_genotype::~em_hmm_genotype()" << endl;
}

void em_hmm_genotype::initialise_model_param()
{
  //m_hap_data.initialise_param(); 
  //m_hap_data.log_param_hap();  
  update_HMM_param();
  log_param(); 
}

void em_hmm_genotype::get_total_stateSpace()
{
 // cout<< "CHK em_hmm_genotype::get_total_stateSpace() START" << endl; 
  gen_combi_states(n_ploidy,n_clusters,m_total_states); 
  num_states = m_total_states.size();

     
  //make the dakota map to facilate gradient computation
    int stcount =0;
    map<int,index_vars> temp1,temp2,temp3;
    for(int icount = 0; icount < n_markers;++icount)
    {
       for(int jcount=0; jcount < n_clusters;++jcount)
       {
	    index_vars temp;
	    temp.marker= icount;
	    temp.cluster = jcount;
	    temp1[stcount]=  temp;
	    ++stcount;
       }
    }
    for(int icount = 0; icount < n_markers;++icount)
    {
       for(int jcount=0; jcount < n_clusters;++jcount)
       {
	    index_vars temp;
	    temp.marker= icount;
	    temp.cluster = jcount;
	    temp2[stcount]=  temp;
	    ++stcount;
       }
    }
     for(int icount = 0; icount < n_markers-1;++icount)
    {
	 index_vars temp;
	 temp.marker= icount;
	 temp.cluster = -1;
	 temp3[stcount]=  temp;
         ++stcount;       
    }
    dakota_map.push_back(temp1);
    dakota_map.push_back(temp2);
    dakota_map.push_back(temp3);
   
    int paramCount =0;
    
 /* for(auto &ip:dakota_map)
  {
     cout <<endl << paramCount << endl;
     for(auto &jp:ip)
     {
        cout << jp.first << ":" << "(" << jp.second.marker << ","<< jp.second.cluster<<") " ;	
    }
    ++paramCount;
  }
   stcount =0;
    for(auto  &it:m_total_states)
    {
      cout << stcount << " " ;
      for(auto jt : it.values)
      {
	  cout << jt << " " ;
      }
      cout << "  Weight:" << it.weight;
      cout << endl;
      ++stcount;
    }
    cout << "State space has " << m_total_states.size() << " states" << endl; */
//cout<< "CHK em_hmm_genotype::get_total_stateSpace() END" << endl; 
}

void em_hmm_genotype::get_current_stateSpace()
{
  //cout<< "CHK em_hmm_genotype::get_current_stateSpace() START" << endl;
  std::string line;
  std::ifstream myReadFile;
 
  myReadFile.open("state_space.txt");
  getline(myReadFile,line); // get line from file
  int count=0,pos;
  
  // Read the statespace from previous instance stored in file state_space.txt
  if(myReadFile.good())
  {
     // myReadFile.open("state_space.txt");
      istringstream buf_stream;
      std::vector<vector<int>> genoDataTemp;
      
      while(myReadFile.good())
      {
	if(line.find("statespace")!=string::npos) // search
	{
	   pos=count;
	}
	getline(myReadFile,line); // get line from file
	++count;
      }
     
      myReadFile.close();
      myReadFile.open("state_space.txt");
      count =0;
      while(myReadFile.good() && count<pos)
      {
	getline(myReadFile,line);
	++count;
      }
      getline(myReadFile,line);
      m_current_states.clear();
      
      for(int icount=0;icount<n_individuals;++icount)
      {
	 getline(myReadFile,line);
	 getline(myReadFile,line);
	 vector<vector<int>> indData;
	 for(int jcount=0;jcount<n_markers;++jcount)
	 {
	   getline(myReadFile,line);	
           buf_stream.str(line);
           genoDataTemp.emplace_back(istream_iterator<int>(buf_stream), istream_iterator<int>());
	   buf_stream.clear();	   
	 }
	 
	 for(int j=0;j<n_markers;++j)
	 {
	    vector<int> temp;
	    for(int i =0; i < num_states;++i)
	      temp.emplace_back(genoDataTemp[j][i]);
	      
	    indData.emplace_back(temp);
	 }
	 genoDataTemp.clear();
	 m_current_states.push_back(indData);
      } 
  }
  
  else // sample the statespace for the very first instance 
  {
   /* std::ofstream myFile;
    int indcnt=0,marker;
    myFile.open("state_space.txt"); 
   
    //we sample separate state space for different individuals
    for(int indCnt=0;indCnt<n_individuals;++indCnt)
    {
      vector<vector<int>> indData;
      for(int markCount=0; markCount < n_markers; markCount++)   
      {
	  vector<int> markerData;	  
	 //sample indices of 10 states according to theirs weights randomly
	  sample(m_total_states,markerData,num_states);	  
	  indData.emplace_back(markerData);
      }
       m_current_states.emplace_back(indData);
    }
    myFile << "statespace 0:" <<  endl<< endl;

    for(auto  &it:m_current_states)
    {
	myFile << "Ind: " << indcnt << endl ;
	for(auto jt : it)
	{   
	    for(auto kt : jt)
	      myFile << kt << " " ;
	  myFile << endl;
	}
	++indcnt;
	myFile << endl;	
     }
 
    myFile.close();*/
  }
  /* 
  //we sample same state space for different individuals
  vector<vector<int>> indData;
   for(int markCount=0; markCount < n_markers; markCount++)   
   {
	vector<int> markerData;	  
	//sample indices of 10 states according to theirs weights randomly
	sample(m_total_states,markerData,num_states);	  
	 indData.emplace_back(markerData);
    }
   for(int indCnt=0;indCnt<n_individuals;++indCnt)
   {
       m_current_states.emplace_back(indData);
   }*/
   
   
     //print values
   /*  cout << "Current state space" << endl;
     int marker =0, indcnt=0;
     for(auto  &it:m_current_states)
     {
        marker =0;
	cout << "Ind: " << indcnt << endl ;
	for(auto jt : it)
	{   
	  cout << "marker: " << marker << " " ;
	    for(auto kt : jt)
	      cout << kt << " " ;
	  ++marker;
	  cout << endl;
	}
	++indcnt;
	cout << endl;
	
     }*/
  //cout<< "CHK em_hmm_genotype::get_current_stateSpace() END" << endl; 
}

void em_hmm_genotype::update_HMM_param()
{	
   //cout<< "CHK em_hmm_genotype::update_HMM_param() START" << endl; 
   clear_variables();  
   m_hap_data.compute_trans_prob_bet_clst();
   compute_jump_prob();  
    
    vector<int> time_count;
    
  for(int marker=0;marker<n_markers;++marker)
    m_clst_given_v .emplace_back(vector<Double>());
	
    for(int indCount=0;indCount<n_individuals;++indCount)
    {
	m_geno_given_clst_v.emplace_back(vector<vector<Double>>());   
	m_clst_given_geno_v.emplace_back(vector<vector<Double>>());
	m_fwd_probs.emplace_back(vector<vector<Double>> ());
	m_bckwd_probs.emplace_back(vector<vector<Double>> ());      
    }
  
    for(int indCount=0;indCount<n_individuals;++indCount)
      for(int marker=0;marker<n_markers;++marker)
      {
	m_fwd_probs[indCount].emplace_back(vector<Double>(num_states));
	m_bckwd_probs[indCount].emplace_back(vector<Double>(num_states));
      }
      
   try
   {
        compute_clust_given_v();  
	compute_geno_given_clst_v(); 
	compute_fwd_prob();  
	compute_bckwd_prob();  
	compute_clst_given_geno_v(); 
	
	/*tbb::tick_count before = tbb::tick_count::now();
	compute_clust_given_v();  
	tbb::tick_count after = tbb::tick_count::now();
	time_count.emplace_back((after - before).seconds() );
		
	before = tbb::tick_count::now();
	compute_geno_given_clst_v(); 
	after = tbb::tick_count::now();
	time_count.emplace_back((after - before).seconds() );
	       
	before = tbb::tick_count::now();
	compute_fwd_prob();  
	after = tbb::tick_count::now();
	time_count.emplace_back((after - before).seconds() );
		
	before = tbb::tick_count::now();
	compute_bckwd_prob();  
	after = tbb::tick_count::now();
	time_count.emplace_back((after - before).seconds() );
		
	before = tbb::tick_count::now();
	compute_clst_given_geno_v(); 
	after = tbb::tick_count::now();
	time_count.emplace_back((after - before).seconds() );*/
	  
   }
   catch(const std::exception& e)
   {
	std::cout << e.what() << '\n';
   }
  // log_param();
 /* std::ofstream myReadFile;
    myReadFile.open("CHKING.txt");
    myReadFile<< "compute_clust_given_v" << " compute_geno_given_clst_v  " << " compute_fwd_prob " 
     " compute_bckwd_prob "<< " compute_clst_given_geno_v "<<  endl ;

     for(auto  &it:time_count)
     {
	   myReadFile<< it << " ";      
     }
   
   //compute_exp_jump_given_geno();
   //compute_jump_given_geno();
    
myReadFile.close(); */

 // cout<< "CHK em_hmm_genotype::update_HMM_param() END" << endl;
}

//p(Jim=j|r)
void em_hmm_genotype::compute_jump_prob()
{
   //cout<< "CHK em_hmm_genotype::compute_jump_prob() START" << endl;
   
   Double jumpProbTemp,jumpProb,dTemp;
  
   int iscaledTimes;
   Double dCopyRecomb; 
   
   for(int jCount=0; jCount < n_markers-1; jCount++)
   {
	vector<Double> markerdata; 
	dTemp= m_hap_data.m_recombinations[jCount]*m_hap_data.m_input.m_physical_distances[jCount];
	//dTemp= 0.0001*m_hap_data.m_input.m_physical_distances[jCount];  
	jumpProbTemp = exp(-1*dTemp);
	
	
	for(int iCount=0;iCount <= n_ploidy;iCount++)
	{
	    jumpProb = binomialProb(n_ploidy,iCount,jumpProbTemp);
	     
	    markerdata.push_back(jumpProb); 
	}
        
        m_jump_prob[jCount] = markerdata;  
    } 
  
  //cout << "em_hmm_genotype::compute_jump_prob() END" << endl;
}

//p(zim|v)
void em_hmm_genotype::compute_clust_given_v()
{
  // cout<< "CHK em_hmm_genotype::compute_clust_given_v() START" << endl;
 
   Double dThresholdvalue = pow(10,-10);
      
  try
   {
        //for initial marker 0 
    
     //for(int indCount=0;indCount<2;++indCount)
     {
          vector<Double>initialstateprobs(num_states);
	   parallel_for(int(0), num_states, [this,&dThresholdvalue,&initialstateprobs](int stateCount)throw()
	  {
	      vector<int> toStVec,tempVec;
	      vector<vector<int>> fromPerms;
	      Double dFinalValue =1.0;
	     // toStVec = m_total_states[stateCount].values;
	       toStVec = m_total_states[stateCount];
	      vector_permutation(toStVec,tempVec,fromPerms);
	      	
	      for(auto  &iter:toStVec)
	      {
		dFinalValue = dFinalValue * m_hap_data.m_alpha[0][iter];
		//dFinalValue = exp(log(dFinalValue)+log(m_hap_data.m_alpha[0][iter]));
	      }
	      dFinalValue = dFinalValue * fromPerms.size();
	      initialstateprobs[stateCount]=(dFinalValue);	 
	      fromPerms.clear();
	    });
	  m_clst_given_v[0]= (initialstateprobs);
   }
     //fromPerms.clear();
  //for markers 1 to n_markers
  
 // for(int indCount=0;indCount<2;++indCount)
  {
	  for(int markCnt=1; markCnt < n_markers; markCnt++)
	  {
	      vector<Double> markerData(num_states);
	      vector<vector<int>> markerTransitions(num_states);
	    //  cout << "markCnt " << markCnt <<endl;
	    
	
	    parallel_for(int(0), num_states, [this,&dThresholdvalue,&markerData,&markerTransitions,&markCnt](int stateCount)throw()
	    {
		  vector<int> toStVec;
		  vector<vector<int>> fromPerms;
		  //toStVec  = m_total_states[stateCount].values;
		   toStVec  = m_total_states[stateCount];
		  Double dTemp,dFinalValue=0.0;
		  vector<int> stateValues;
		  vector<int> frmStVec,tempVec;
		
		  for(int fromCount=0;fromCount< num_states;++fromCount)
		  {
		      // cout <<fromCount << " ";
			//frmStVec = m_total_states[fromCount].values;	 
		        frmStVec = m_total_states[fromCount];
			//get the permutations of the from state vector
			vector_permutation(frmStVec,tempVec,fromPerms);
			//dFinalValue = dFinalValue + m_hap_data.compute_trans_prob_bet_clst_tuples(markCnt-1,fromPerms,toStVec);
			dTemp =  m_hap_data.compute_trans_prob_bet_clst_tuples(markCnt-1,fromPerms,toStVec);
			fromPerms.clear();
			dFinalValue += dTemp;
			//if(dTemp >dThresholdvalue)
			{
			   stateValues.emplace_back(fromCount);
			}
		}
		    
		    markerData[stateCount]=dFinalValue;
		    //markerTransitions[stateCount]=stateValues;
	      });
	        m_clst_given_v[markCnt]=markerData;   
		//m_current_transitions.emplace_back(markerTransitions);
	    }
    }
 }
  catch(const std::exception& e)
  {
	std::cout << e.what() << '\n';
  }
   
 // cout<< "CHK em_hmm_genotype::compute_clust_given_v() END" << endl; 
}

 
void em_hmm_genotype::compute_geno_given_clst_v()
{
 // cout<< "CHK em_hmm_genotype::compute_geno_given_clst_v() START" << endl;

  vector<Double> thetaCurMarker, clustCurState;
  vector<int> haplo_vector,temp;//currentMark_states
  vector<vector<int>> genoCurInd, permuted_haplotypes;
  
  vector< vector<std::pair<int,int>>> Map_Haplos_states;
  vector<Double>clustCurMarker;

  try
  {
         for(int indCount=0;indCount<n_individuals;++indCount)
	 {
	      genoCurInd = m_hap_data.m_input.m_genotypes[indCount];   
	
	      for(int markCount=0; markCount < n_markers; ++markCount)  
	      { 
		    vector<Double> markerData(num_states);
		    thetaCurMarker = m_hap_data.m_theta[markCount];
		      clustCurMarker = m_clst_given_v[markCount];
		    //get the haplo vector.For eg if geno is 0 for tetraploids,  haplo_vector= {0, 0, 0, 0} as given in the inout file
		    haplo_vector = genoCurInd[markCount];  
		    vector_permutation(haplo_vector,temp,permuted_haplotypes);
		   
		    //for(int stateCount=0; stateCount < num_states;++ stateCount)
		    parallel_for(int(0), num_states, [this,&permuted_haplotypes,&markerData,&thetaCurMarker,&clustCurMarker](int stateCount)throw()
		    {
			  long int ihaplo;
			  Double final_value =0.0,dTemphaplo,perm_values,dThetaOne; 
			  vector<int>  cur_state;
			  vector< vector<std::pair<int,int>>> Map_Haplos_states;
			  // cur_state = m_total_states[currentMark_states[stateCount]].values;
			 // cur_state = m_total_states[stateCount].values;
			  cur_state = m_total_states[stateCount];
			  // map the state and haplo combi in unique way, for eg haplo={0, 0, 1, 1} state={2,3,3,4} map= {(0, 2),(0, 3),(1, 3),(1, 4)}
			  get_map_vectors(permuted_haplotypes, cur_state, Map_Haplos_states);
			
			  for(auto &kvp:Map_Haplos_states)
			  {
			      perm_values=1.0;
			    
			      for (auto clst_it:kvp)
			      {
				  dThetaOne = thetaCurMarker[clst_it.second];
				  ihaplo =   clst_it.first;
			      
				// Equation (1) in Scheet & Stephans 2006  article
				//dTemphaplo = pow(dThetaOne,ihaplo)* pow((1.0-dThetaOne),(1.0-ihaplo));
				dTemphaplo = pow(dThetaOne,ihaplo);
				dThetaOne = 1.0-dThetaOne ; 
				ihaplo = 1.0-ihaplo;
				dTemphaplo *= pow(dThetaOne,ihaplo);		    
				perm_values = perm_values*dTemphaplo;		 
			    }
			    //sum of all possible permutations     
			    final_value +=  perm_values ;	    
			}
			Map_Haplos_states.clear();
		      // final_value = exp(log(final_value)+log(clustCurState[stateCount]));
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

     
void  em_hmm_genotype::util_hmm_fwd(int marker,int state,int njumps,vector<vector<int>> &jumpData, 
			  vector<Double>&alpha_data)
{
 // cout<< "CHK em_hmm_genotype::util_hmm_fwd() START" << endl;
  vector<Double> AlphaMarker;
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

   //data from PARTICLE FILTER
    //currentMark_trans = m_current_transitions[marker-1][state];
 
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
	       //REMOVE those states filtered out by PARTICLE FILTER
		  if(std::find(currentMark_trans.begin(),currentMark_trans.end(),state_jump)!= currentMark_trans.end())
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
void em_hmm_genotype::util_hmm_bckwd(int marker,int state,int njumps,vector<vector<int>> &jumpData,
			vector<vector<Double>>&alpha_data)
{
  //cout<< "CHK em_hmm_genotype::util_hmm_bckwd() START" << endl;
  
  vector<Double> AlphaMarker;
  int state_jump;
  Double dalpha_temp;
  vector<int> input_comb,copy_Curstate,currentMark_trans;
   
    vector< vector<int>> ::iterator itTemp;
   
   // WE FIRST work to obtain the jumpData
   vector<int> curState= m_total_states[state];
   vector<int> gotTemp;  
   vector< vector<int>> T,Tfinal;
   vector<vector< vector<int>>> Talpha;
     // fetch all the possible combinations of clusters with given number of jumps in current state
   vector_combination(curState,gotTemp, 0, njumps,T);

   // fetch data provided by PARTICLE FILTER
     //currentMark_trans = m_current_transitions[marker][state];
   
   for(int i=0;i<n_clusters;++i)
      input_comb.emplace_back(i);
   vector<int>::iterator pos;
   
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
	    itTemp = std::find( m_total_states.begin(), m_total_states.end(), iter);
	    
	    if(itTemp!=m_total_states.end())
	    {
	         state_jump = std::distance( m_total_states.begin(),itTemp);
	       //REMOVE those states filtered out by PARTICLE FILTER
		 // if(std::find(currentMark_trans.begin(),currentMark_trans.end(),state_jump)!= currentMark_trans.end())
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
 // cout<< "CHK em_hmm_genotype::compute_fwd_prob()" << endl;  
  vector<Double> genoGivenClst_cur; 
  vector<vector<Double>>  genoCurInd;
  Double dTempFwdProb; 
  vector<Double> AlphaCurMarker;
    
  AlphaCurMarker = m_hap_data.m_alpha[0];
  for(int indCount=0;indCount<n_individuals;++indCount)
   {
	// indStates = m_current_states[indCount];
	genoGivenClst_cur = m_geno_given_clst_v[indCount][0];
	parallel_for(int(0), num_states, [this,&genoGivenClst_cur,&AlphaCurMarker,&indCount](int stateCount)throw()    
	{
	      vector<int> curState= m_total_states[stateCount];
	      Double dTempFwdProb =1.0;
	      //for(vector<int>::iterator j = curState.values.begin();j != curState.values.end();j++)
	      for(vector<int>::iterator j = curState.begin();j != curState.end();++j)
	      {
		    dTempFwdProb = dTempFwdProb *AlphaCurMarker[(*j)];
	      }
	      dTempFwdProb = dTempFwdProb * genoGivenClst_cur[stateCount];
	      m_fwd_probs[indCount][0][stateCount]=(dTempFwdProb);
	  });
  }
  
  vector<Double> jumpProbPrevMarker;
  vector<vector<Double>> FwdProbCurInd;//,alphaValTostates;

for(int markCnt=1; markCnt < n_markers ;markCnt++)
{
      jumpProbPrevMarker = m_jump_prob[markCnt-1];
      
      parallel_for(int(0), num_states, [this,&jumpProbPrevMarker,&AlphaCurMarker,
			    &markCnt](int toStCnt)throw()    
      {
		     
	  vector<int> toStVec,currentMark_trans;
	  vector<Double>FwdProbPrevMarker;
	  Double jumpVal,dsumalphaVal,dTemp,dFwdProbVal = 0.0;		    
	  		      
	  toStVec = m_total_states[toStCnt];
	  //currentMark_trans = m_current_transitions[markCnt-1][toStCnt];	
		 		     
	  vector<vector<vector<int>>> JumpDatavec;
	  vector<vector<Double>> alphaValTostates;
	  vector<vector<int>> states_jump;
	  vector<Double> alpha_jump;
	  //required for case where njumps=n_ploidy
	  Double dalphaValstate=1.0;	  
	  for(int ncount=0;ncount<n_ploidy;++ncount)
	  {
		 dalphaValstate *= AlphaCurMarker[toStVec[ncount]];
           }
	    //for jumps from 1 to n_ploidy, get jumpProbability 
	   //and corresponding alphaval contributionbased on num jumps
	  for(int jumpCount=1; jumpCount < n_ploidy;++jumpCount)
	  {
		vector<vector<int>> curJumpData;
		vector<Double>alphaValTostates_Temp;
		util_hmm_fwd(markCnt,toStCnt,jumpCount,curJumpData,alphaValTostates_Temp);			   
		JumpDatavec.emplace_back(curJumpData);	
		alphaValTostates.emplace_back(alphaValTostates_Temp);
	  }
	
	for(int indCount=0;indCount<n_individuals;++indCount)
	  {
		FwdProbPrevMarker = m_fwd_probs[indCount][markCnt-1];			    
		             
		//num of jumps is zero
		/*auto position =std::find(currentMark_trans.begin(), currentMark_trans.end(), toStCnt);
		if(position != currentMark_trans.end())
		    dFwdProbVal = FwdProbPrevMarker[std::distance(currentMark_trans.begin(),position)] * jumpProbPrevMarker[0];*/
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
		  //for(int ncount=0;ncount<currentMark_trans.size();++ncount)
		  for(int ncount=0;ncount<num_states;++ncount)
		  {
			 dTemp+= FwdProbPrevMarker[ncount];
		  }
		  dFwdProbVal += dTemp *jumpProbPrevMarker[n_ploidy]*dalphaValstate;		 
		 dFwdProbVal *= m_geno_given_clst_v[indCount][markCnt][toStCnt];
		 m_fwd_probs[indCount][markCnt][toStCnt]=dFwdProbVal;		
		/* alphaValTostates.clear();
		 JumpDatavec.clear();*/ 
	    }
    });
 }
 //  cout << "em_hmm_genotype::compute_fwd_prob() END" << endl;
}
void em_hmm_genotype::compute_bckwd_prob()
{
  // cout << "em_hmm_genotype::compute_bckwd_prob() START" << endl;
 
     //The backward prob at the Final marker is 1 
   for(int indCount=0;indCount<n_individuals;++indCount)
   {
	 for(int i=0;i < num_states;++i)
	 {
	      m_bckwd_probs[indCount][n_markers-1][i]=(1.0);
	  }
    }
    
   vector<Double> BckwrdProbNxtMarker;
    vector<Double> jumpProbNxtMarker;
     vector<Double> AlphaNXtMarker;
    for(int markCnt = n_markers-2 ; markCnt >= 0 ;--markCnt)
    {
	jumpProbNxtMarker = m_jump_prob[markCnt];
	AlphaNXtMarker = m_hap_data.m_alpha[markCnt+1];
	
	parallel_for(int(0), num_states, [this,&jumpProbNxtMarker,&AlphaNXtMarker,
				&markCnt](int toStCnt)throw()    
	{
	      vector<vector<int>> states_jumpTemp;
	      vector<vector<Double>> alpha_jumpTemp;
	      vector<int> toStVec,currentMark_trans,curJumpstates;
	      vector<Double>BckwrdProbNxtMarker;
	      vector<Double>geno_given_clst_v;
	      Double jumpVal,dsumalphaVal,dTemp,dBwdProbVal = 0.0;		    
	      		      
	     // toStVec = m_total_states[toStCnt].values;
	      toStVec = m_total_states[toStCnt];
	      //currentMark_trans = m_current_transitions[markCnt][toStCnt];	
		 		     
	      vector<vector<vector<int>>> JumpDatavec;
	      vector<vector<vector<Double>>> alphaValTostates;
	     
	       //for jumps from 1 to n_ploidy, get jumpProbability 
	   //and corresponding alphaval contributionbased on num jumps
	  for(int jumpCount=1; jumpCount < n_ploidy;++jumpCount)
	  {
		vector<vector<int>> curJumpData;
		vector<vector<Double>>alphaValTostates_Temp;
		util_hmm_bckwd(markCnt,toStCnt,jumpCount,curJumpData,alphaValTostates_Temp);			   
		JumpDatavec.emplace_back(curJumpData);	
		alphaValTostates.emplace_back(alphaValTostates_Temp);
	  }
	for(int indCount=0;indCount<n_individuals;++indCount)
	  {
	        dBwdProbVal=0.0;
		BckwrdProbNxtMarker = m_bckwd_probs[indCount][markCnt+1];
		geno_given_clst_v = m_geno_given_clst_v[indCount][markCnt+1];
		//num of jumps is zero
		/*auto position =std::find(currentMark_trans.begin(), currentMark_trans.end(), toStCnt);
		if(position != currentMark_trans.end())
		{
		    dBwdProbVal = BckwrdProbNxtMarker[std::distance(std::begin(currentMark_trans),position)]
		       *jumpProbNxtMarker[0]* geno_given_clst_v[toStCnt];
		}*/
		dBwdProbVal = BckwrdProbNxtMarker[toStCnt]
		       *jumpProbNxtMarker[0]* geno_given_clst_v[toStCnt];
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
		// for(int ncount=0;ncount < currentMark_trans.size();++ncount)
		for(int ncount=0;ncount < num_states;++ncount)
		{
		      dTemp = 1.0;
		      for(int icount = 0;icount<n_ploidy;++icount)
			  dTemp *=AlphaNXtMarker[icount];
		      dsumalphaVal += (dTemp *geno_given_clst_v[ncount]*BckwrdProbNxtMarker[ncount]);			
		}
		dsumalphaVal *=  jumpProbNxtMarker[n_ploidy];
		dBwdProbVal +=  dsumalphaVal;
		m_bckwd_probs[indCount][markCnt][toStCnt]=dBwdProbVal;	  
	  }
	  
	});
    }
      
 // cout << " em_hmm_genotype::compute_bckwd_prob() END " << endl;
}

Double em_hmm_genotype::compute_geno_given_v(int ind, int mark)
{
   //cout << "em_hmm_genotype::compute_geno_given_v() Start" << endl;
  
    Double dTempmark=0.0;
    int stateIndex;
    vector<Double > genoGivenClst_cur,fwdProb;
 
    int iscaledTimes = 0;
    
    if(mark >=n_markers)
    {
	return(-1.0);
    }
    
    fwdProb = m_fwd_probs[ind][mark];
    
    for(std::vector<Double>::iterator i=fwdProb.begin();i!=fwdProb.end();i++) 
    {
	dTempmark = dTempmark + ((*i));
    }
   
   //cout << "em_hmm_genotype::compute_geno_given_v()  End" << endl; 
    
    return(dTempmark); 
}

void em_hmm_genotype::compute_clst_given_geno_v()
{
    // cout<< "CHK compute_clst_given_geno_v Start" << endl;  
      vector<Double> FwdProbCurMark, BckwrdProbCurMark;  
      vector<int> currentMark_states,stateVec,tempVec;
      vector<vector<int>> perms_frmState,indData;
      int curState,numPerms;
      Double dTemp;
      
      for(int indCount = 0; indCount < n_individuals; ++indCount)
      {
	  //vector<vector<Double>>  individualData;  
	// indData = m_current_states[indCnt];
	  for(int markCnt=0; markCnt < n_markers ;++markCnt)
	  {
		  vector<Double> markerData;
		  FwdProbCurMark = m_fwd_probs[indCount][markCnt];
		  BckwrdProbCurMark = m_bckwd_probs[indCount][markCnt];
		// currentMark_states = indData[markCnt];
	    
		  for(int stateCount=0; stateCount < num_states; ++stateCount)
		  {
		      //curState = currentMark_states[stateCount] ;
		    // stateVec = m_total_states[stateCount].values;
		  
		      dTemp = FwdProbCurMark[stateCount]*BckwrdProbCurMark[stateCount];
		      markerData.emplace_back(dTemp);
		  }
		  
		  m_clst_given_geno_v[indCount].push_back(markerData);
	    }
	  
	    //m_clst_given_geno_v.push_back(individualData);
	  //m_clst_given_geno_v.insert(m_clst_given_geno_v.begin()+indCnt,individualData);
 }
  
    /*  std::vector< vector <Double>> IndData = m_clst_given_geno_v[0];  
      
      for(std::vector< vector<Double>>::iterator i = IndData.begin(); i != IndData.begin()+2;i++)
      {
	    for(std::vector<Double>::iterator j= (*i).begin(); j!=(*i).end() ;j++) 
	    {
		cout <<  (*j) << "  " ;
	    }
	  
	    cout << endl << " ************** **************** ************* ************ *********** " << endl;
       }*/
    
 // cout<< "CHK compute_clst_given_geno_v End" << endl;
    
}
  
void em_hmm_genotype::log_param()
{
    std::ofstream myReadFile;
    int marker =0, indcnt=0,state;    
    myReadFile.open("CHKING2.txt");
    
     myReadFile << endl<<"***** compute_clust_given_v*****"  << endl;
  
  for(auto  &it:m_clst_given_v)
  { 
	myReadFile<<endl << "marker: " << marker << endl ;
	state=0;
	for(auto jt : it)
	{   
	        myReadFile <<  jt<< " " ;
	 }
	 ++marker;
  } 
 myReadFile << endl <<" ********compute_geno_given_clst_v*********"  << endl;
indcnt=0;
for(auto  &it:m_geno_given_clst_v)
  {
    myReadFile << "IND :" << indcnt << endl;
    marker =0;
      for(auto jt : it)
      {
	myReadFile << "Marker "<< marker <<endl;
	for(auto kt : jt)
	{
	  myReadFile << kt << " ";
	}
	++marker;
	myReadFile << endl;
      }
     ++indcnt;
  }
 myReadFile <<endl << " **********Fwd probs**********"  << endl;
      indcnt=0;  
      for(auto &kvp: m_fwd_probs)
      {
	 myReadFile << "Ind " << indcnt << endl ;
	 marker = 0;
	  for(auto &kv:kvp) 
	  {
	      myReadFile << "Marker " << marker <<endl;
	      for(auto p:kv)
		myReadFile << p <<  " ";  
	      ++marker;
	      myReadFile << endl;
	  }
	 ++indcnt;
       }
       
        myReadFile <<endl << " **********Bckwd probs**********"  << endl;
      indcnt=0;  
      for(auto &kvp: m_bckwd_probs)
      {
	 myReadFile << "Ind " << indcnt << endl ;
	 marker = 0;
	  for(auto &kv:kvp) 
	  {
	      myReadFile << "Marker " << marker <<endl;
	      for(auto p:kv)
	       myReadFile << p <<  " ";  
	      ++marker;
	      myReadFile << endl;
	  }
	 ++indcnt;
       }
    
       myReadFile <<endl << " **********compute_clst_given_geno_v**********"  << endl;
      indcnt=0;  
      for(auto &kvp: m_clst_given_geno_v)
      {
	 myReadFile << "Ind " << indcnt << endl ;
	 marker = 0;
	  for(auto &kv:kvp) 
	  {
	      myReadFile << "Marker " << marker <<endl;
	      for(auto p:kv)
	       myReadFile << p <<  " ";  	      
	      ++marker;
	      myReadFile << endl;
	  }
	 ++indcnt;
       }
 /*myReadFile << endl <<" ********chosen transitions*********"  << endl;
indcnt=0;
for(auto  &it:m_current_transitions)
{
    myReadFile << "Marker "<< marker << " ";
    state =0;
    
    for(auto jt : it)
    {
	myReadFile << "State "<< state << " ";
	for(auto kt : jt)
	{
	    myReadFile << kt << " ";
	}
	++state;
	myReadFile << endl;
      }
     ++marker;
}*/

 myReadFile.close();
}

void em_hmm_genotype::clear_variables()
{
  //cout<< "CHK em_hmm_genotype::clear_variables() START" << endl;
  m_hap_data.m_trans_prob_bet_clst.clear();
  m_jump_prob.clear();
  m_geno_given_clst_v.clear();
  m_fwd_probs.clear();
  m_bckwd_probs.clear();
  m_clst_given_geno_v.clear();
  m_clst_given_v.clear();
  //m_current_transitions.clear();
   //cout<< "CHK em_hmm_genotype::clear_variables() END" << endl;
}

void em_hmm_genotype::update_statespace()
{
   vector<vector<vector<int>>> stateSpaceUpdated;  
  
  /* m_part_filter.filter_run(m_fwd_probs,m_bckwd_probs,
			      m_current_states,m_total_states,stateSpaceUpdated);*/
   m_current_states.clear();	      
   m_current_states = stateSpaceUpdated;
}

//p(zim|gim,jim v*}*{log[p(gim|zim, jim, v)]+log[p(zim|v,jim)]}
Double em_hmm_genotype::func_eval_local()
{
    cout<< "CHK em_hmm_genotype::func_eval_local() START" << endl; 
    update_HMM_param();
    Double functionEval = 1.0, dIndValue,dMarkValue,dClstGgenoV,dGenoGclustV,dClustGv; 
    vector<Double> Geno_given_Clst_V_CurMrk,Clst_given_Geno_V_CurMrk,Clst_given_V;
    std::ifstream myReadFile;
    int dScaleFactor=0;
  //  Double dscaler=pow(10,30);
   // Double dThreshold=pow(10,-40);
    int scaler =0;
    string Strbuffer;
   
    for(int indCnt = 0; indCnt< n_individuals;indCnt++)  
    {
	dIndValue = 1.0;
	
	for(int iCurMarker = 0; iCurMarker < n_markers; iCurMarker++)
	{
	     dMarkValue=0.0;
	     
	     Clst_given_Geno_V_CurMrk = m_clst_given_geno_v[indCnt][iCurMarker];//p(zim|gim, v*)
	     Geno_given_Clst_V_CurMrk = m_geno_given_clst_v[indCnt][iCurMarker];//p(gim|zim, jim, v)
	     //Clst_given_V = m_clst_given_v[iCurMarker]; // p(zim|jim,v)
	     
	     for(int iStateCount=0; iStateCount < num_states; iStateCount++)	
	     {
	           dClstGgenoV = Clst_given_Geno_V_CurMrk[iStateCount];       
		   dGenoGclustV = Geno_given_Clst_V_CurMrk[iStateCount];	
		   //dClustGv = Clst_given_V[iStateCount];
		   dMarkValue = dMarkValue + dClstGgenoV*(log(dGenoGclustV));//+log(dClustGv));
	     }
	    dIndValue = dIndValue*dMarkValue;
	}
	//functionEval = functionEval + dIndValue;
	
	functionEval *=  dIndValue;
    }
    cout<< "CHK em_hmm_genotype::func_eval_local() B4: " << functionEval<< endl << dScaleFactor << endl; 
    functionEval = functionEval*-1.0;
   
   /*myReadFile.open("Phase_param.txt");
   
  while(myReadFile)
   {
      getline(myReadFile,line);
	
      if(line.find("statespace") != string::npos) // search
	++count;
   } 
    myReadFile.close();
     update_statespace();
   //m_params_file << functionEval<< endl;
  // m_params_file << ".............................................................................." << endl;
   //m_params_file << "statespace " << count << ":" << endl << endl ; 
   log_param();
   */
   ifstream scale_file("scaling.in");
  
   if(scale_file.is_open())
   {
     getline(scale_file,Strbuffer);
     scaler = std::stoi(Strbuffer) ;
    
     scale_file.close();
     
     if(scaler>0)
     {
	  while(scaler>0)// && abs(functionEval)<dThreshold)
	  {
		functionEval = functionEval*scaling_factor;
		--scaler;	    
	  }
      }
      cout << "Scaler " <<  scaler << endl; 
       cout<< "CHK em_hmm_genotype::func_eval_local() Aft: " << functionEval<< endl; 
  }
  
    return functionEval;
}

//p(zim|gim,jim v*}*{1/p(gim|zim, jim, v)}{d/dtheta_km}[p(gim|zim,Jim,theta_km)p(zim|alpha_km,rm)]
Double em_hmm_genotype::get_theta_gradient(int iCluster, int iCurMarker)
{
    //cout<< "CHK em_hmm_genotype::get_theta_gradient() START" << endl; 
    Double df_theta_km,dtot_IndValue ;
    int iCurCluster = iCluster;
    int loopEnd;
    df_theta_km = 0.0;
    vector<Double> Clst_given_Geno_V_CurMrk, Geno_given_Clst_V_CurMrk,Clst_given_V ,thetaCurMarker;
    vector<int> haplo_vector,states_relevant,tempVec,toStVec;
    vector<vector<int>> permuted_haplotypes;    
    
    combinable<Double> gradtheta_state([]() { return 0.0; });     
    Clst_given_V = m_clst_given_v[iCurMarker]; // p(zim,jim|v)
    thetaCurMarker = m_hap_data.m_theta[iCurMarker];
    
    
    //get all the states with our cluster  'iCurCluster' at the current marker iCurMarker ; rest of them become zero during differentiation.        
    for(int iStateCount=0; iStateCount < num_states; ++iStateCount)
    {
         toStVec = m_total_states[iStateCount];
	 if(std::find(toStVec.begin(), toStVec.end(), iCurCluster) != toStVec.end())
	 {
	      states_relevant.push_back(iStateCount);
	  }
     }
     loopEnd = states_relevant.size();
     

   for(int indCnt = 0; indCnt< n_individuals;indCnt++)
   {
	 Clst_given_Geno_V_CurMrk = m_clst_given_geno_v[indCnt][iCurMarker];//p(zim|gim, Jim,v*)
	 Geno_given_Clst_V_CurMrk =  m_geno_given_clst_v[indCnt][iCurMarker];//p(gim|zim, Jim,v)
	 haplo_vector = m_hap_data.m_input.m_genotypes[indCnt][iCurMarker];
	 dtot_IndValue = 0.0; 	
	
	 //we permute the haplotypes,  but not the states
           vector_permutation(haplo_vector,tempVec,permuted_haplotypes);
	      // compute {d/dtheta_km}[p(gim|zim, jim, v)], which is {d/dtheta_km}[log (p(gim|zim,theta_km)*p(zim|alpha_km,Jim,rm))]
	   parallel_for(int(0), loopEnd, [this,&Clst_given_Geno_V_CurMrk,&Geno_given_Clst_V_CurMrk,&states_relevant, 
	      &thetaCurMarker, &permuted_haplotypes,&Clst_given_V,&gradtheta_state,&iCurCluster](int stateCount)throw()    
	    {
	        vector<int> cur_state;
		vector<Double> variable_vector_var,variable_vector_const;
		int icounter,ihaplo;
		vector< vector<std::pair<int,int>>> Map_Haplos_states;
		Double const_value,var_value,dtot_stateValue=0.0,dTemphaplo,dThetaOne, dClstGgenoV,dGenoGclustV,dTempvar; 
		
	        dClstGgenoV = Clst_given_Geno_V_CurMrk[states_relevant[stateCount]];       
		dGenoGclustV = Geno_given_Clst_V_CurMrk[states_relevant[stateCount]];					
       	        cur_state = m_total_states[states_relevant[stateCount]];
		
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
		      dtot_stateValue = dtot_stateValue + dTempvar;			       
		      		     
		}
	       // include the constant part of the entire expression
		  dtot_stateValue  = dtot_stateValue * (dClstGgenoV/dGenoGclustV); 	
		   // sum of all derivative for all the relevant states
		  gradtheta_state.local() += dtot_stateValue;
		 // dtot_IndValue = dtot_IndValue + dtot_stateValue;   
		 // permuted_haplotypes.clear();
		 // Map_Haplos_states.clear();
	    });
	//df_theta_km +=dtot_IndValue;      
      df_theta_km +=gradtheta_state.combine(plus<Double>());
      gradtheta_state.clear();
 }
 
  //cout<< "CHK em_hmm_genotype::get_theta_gradient() END" << endl; 
  return df_theta_km;
}

Double  em_hmm_genotype::diff_eval_local(int index)
{
  //cout<< "CHK em_hmm_genotype::diff_eval_local() START" << endl; 
 
  int iCurMarker,iCurcluster;
  Double  grad_value = -1.0;

  if(index < (n_markers*n_clusters))    
  {
    /*if(index == 0 )
    {
      cout << endl<<" theta grads" << endl;
    }*/
    iCurMarker = dakota_map[0][index].marker;
    iCurcluster = dakota_map[0][index].cluster;
    cout << index << " " << iCurMarker << " " << iCurcluster << endl;
    grad_value = get_theta_gradient(iCurcluster,iCurMarker);
  }

  else if((index < (2*n_markers*n_clusters)) && (index >=(n_markers*n_clusters)))
  {
    /* if(index == (n_markers*n_clusters) )
    {
      cout << endl<<" alpha grads" << endl;
    }*/
    iCurMarker = dakota_map[1][index].marker;
    iCurcluster = dakota_map[1][index].cluster;
    grad_value = get_alpha_gradient(iCurcluster,iCurMarker);
  }

  else
  {
    /*if(index == (n_markers*n_clusters*2) )
    {
      cout << endl<<" Recomb grads" << endl;
    }*/
    iCurMarker = dakota_map[2][index].marker;
    grad_value = get_recomb_gradient(iCurMarker+1);
  }
  
 return grad_value;
// cout<< "CHK em_hmm_genotype::diff_eval_local() END" << endl; 
}

//differentiating wrt to alpha_km
// 2*p(zim|gim,jim,v*)*{1/p(zim|Jim,alpha_km,rm)}{d/dalpha_km}[p(zim|Jim,alpha_km,rm)]
Double  em_hmm_genotype::get_alpha_gradient(int iCluster,int iCurMarker)
{
 //   cout<< "CHK em_hmm_genotype::get_alpha_gradient() START" << endl; 
  //tbb::tick_count before = tbb::tick_count::now();
  Double df_alpha_km=0.0,dtot_IndValue;
  int iCurCluster = iCluster,loopEnd;
  vector<int> toStVec,states_relevant;
  vector<Double> Clst_given_Geno_V_CurMrk,Clst_given_V;
  combinable<Double> gradalpha_state([]() { return 0.0; });   
  //get all the states with our cluster  'iCurCluster' at the current marker iCurMarker ; rest of them become zero during differentiation.        
  for(int iStateCount=0; iStateCount < num_states; ++iStateCount)
  {
	toStVec = m_total_states[iStateCount];
	if(std::find(toStVec.begin(), toStVec.end(), iCurCluster) != toStVec.end())
	{
	    states_relevant.emplace_back(iStateCount);
	}
   }
   loopEnd = states_relevant.size();
    Clst_given_V =m_clst_given_v[iCurMarker];
   if (iCurMarker > 0) //differentiate p(zim|v) at m >0,  which is function for Transition probs as Eq 9 in Scheet & Stephans 2006
  {
        //vector<Double> transProb_curMarker =   m_hap_data.m_trans_prob_bet_clst[iCurMarker-1];
	for(int indCnt = 0; indCnt< n_individuals;++indCnt)  
        {
	   Clst_given_Geno_V_CurMrk = m_clst_given_geno_v[indCnt][iCurMarker];
	    
	    parallel_for(int(0), loopEnd, [this,&iCurCluster,&states_relevant,&iCurMarker,&gradalpha_state,
	      &Clst_given_V,&Clst_given_Geno_V_CurMrk ](int stateCount)throw()
	    {
	          vector<int> cur_stateVec,frm_stateVec,tempVec;
		  vector<vector<int>> permuted_states;
		  vector< vector<std::pair<int,int>>> Map_fromTo_states; 
		  int icounter;
		  Double const_value, var_value,dTemp,dcurStateValue =0.0,dTempvar,dTempk;
		  vector<Double> variable_vector_const,variable_vector_var;
		  cur_stateVec = m_total_states[states_relevant[stateCount]];
		  
		   for(int fromCount=0; fromCount< num_states;fromCount++)
		   {
			  frm_stateVec = m_total_states[fromCount];
		          vector_permutation(frm_stateVec,tempVec,permuted_states);   
		          get_map_vectors(permuted_states, cur_stateVec, Map_fromTo_states);
			  
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
			    dcurStateValue = dcurStateValue+dTemp;	
		      }
		       Map_fromTo_states.clear();
		       permuted_states.clear();
		  }
		  
		  // include constatnt term of the entire expression 
	      dTempvar = 2*Clst_given_Geno_V_CurMrk[states_relevant[stateCount]];
	      dTempk = Clst_given_V[states_relevant[stateCount]];
	      const_value = (dTempvar/dTempk);
	
	      // const_value = const_value * n_individuals;
	       dcurStateValue = dcurStateValue*const_value;
	       //dtot_IndValue = dtot_IndValue + dcurStateValue; 
	        gradalpha_state.local() += dcurStateValue;
	    });
	    //df_alpha_km = df_alpha_km + dtot_IndValue; 
	      df_alpha_km +=gradalpha_state.combine(plus<Double>());
	      gradalpha_state.clear();
	}
  }
  else // marker =0  differentiate p(zim|v) at m=0,  which is function for Initial state probs as Eq 9 in Scheet & Stephans 2006
  {
	for(int indCnt = 0; indCnt< n_individuals;++indCnt)  
        {
	      parallel_for(int(0), loopEnd, [this,&iCurCluster,&states_relevant,&iCurMarker,&gradalpha_state,&Clst_given_V](int stateCount)throw()
	      {
		  vector<int> cur_stateVec;
		  cur_stateVec = m_total_states[states_relevant[stateCount]];
		  Double const_value =1.0, var_value =0.0,dTemp;
		  //vector_permutation(cur_stateVec,tempVec,permuted_states);
		
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
		  //dtot_IndValue = dtot_IndValue + dTemp;
		  gradalpha_state.local() += dTemp;
		  //permuted_states.clear();   
	     });
	      // states_relevant.clear();
	      df_alpha_km += gradalpha_state.combine(plus<Double>());
	      gradalpha_state.clear();
	}
  }
  states_relevant.clear();
  //cout<< "CHK em_hmm_genotype::get_alpha_gradient() END" << endl; 
  //tbb::tick_count after = tbb::tick_count::now();
//cout << (after - before).seconds()<< endl;
  return(df_alpha_km);
}

//differentiating wrt to r_m,
// p(zim|gim,jim, v*)*2*{1/p(zim|jim,alpha_km,rm)}{d/dr_m}[p(zim|jim,alpha_km,rm))]
Double  em_hmm_genotype::get_recomb_gradient(int iCurMarker)
{
   // cout << "em_hmm_genotype::get_recomb_gradient() START" << endl;  
    Double df_r_m=0.0;
    vector<Double> alphaCurMarker,Clst_given_V,Clst_given_Geno_V_CurMrk;
    combinable<Double> gradrecomb_state([]() { return 0.0; });   
    if (iCurMarker > 0) 
    {
	  alphaCurMarker = m_hap_data.m_alpha[iCurMarker];
	  Clst_given_V =m_clst_given_v[iCurMarker];
	  for(int indCnt = 0; indCnt< n_individuals;++indCnt)  
          {
	        Clst_given_Geno_V_CurMrk = m_clst_given_geno_v[indCnt][iCurMarker];
		
		parallel_for(int(0), num_states, [this,&iCurMarker,&alphaCurMarker,&Clst_given_Geno_V_CurMrk,
		&gradrecomb_state,&Clst_given_V](int iStateCount)throw()
		{
		      vector<int> cur_stateVec,frm_stateVec,tempVec;
		      Double dcurStateValue =0.0,var_value,dTemp,dTempvar,dTempk,const_value;		     
		      vector<vector<int>> permuted_states;
		      vector< vector<std::pair<int,int>>> Map_fromTo_states;
		      vector<Double> variable_vector_const,variable_vector_var;
		      int icounter;
		      cur_stateVec = m_total_states[iStateCount];
		      for(int fromCount=0; fromCount< num_states;fromCount++)
		      {
			    frm_stateVec = m_total_states[fromCount];
			    vector_permutation(frm_stateVec,tempVec,permuted_states);   
			    get_map_vectors(permuted_states, cur_stateVec, Map_fromTo_states);
			   
			    //differentiating wrt to r_m,  the transition probs as Eq 9 in Scheet & Stephans 2006
			    for(auto &kvp:Map_fromTo_states)
			    {
				    var_value =0.0;
				    for(auto from_to:kvp)
				    {
					  dTemp = exp(-1*m_hap_data.m_recombinations[iCurMarker-1]*m_hap_data.m_input.m_physical_distances[iCurMarker-1]);
					  //dTemp = exp(-1*0.0001*m_hap_data.m_input.m_physical_distances[iCurMarker-1]);
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
				    dcurStateValue = dcurStateValue+var_value;
			    }
			    Map_fromTo_states.clear();
			    permuted_states.clear();
		      }
			//include constatnt term of the entire expression		
		      dTempvar = 2*Clst_given_Geno_V_CurMrk[iStateCount];
		      dTempk = Clst_given_V[iStateCount];
		      const_value = (dTempvar/dTempk);		
		      //const_value = const_value * n_individuals;
		      dcurStateValue = dcurStateValue*const_value;
		      //dtot_IndValue = dtot_IndValue + dcurStateValue; 
		      gradrecomb_state.local() += dcurStateValue;
	  });
	 // df_r_m  += dtot_IndValue;  
	df_r_m += gradrecomb_state.combine(plus<Double>());
	gradrecomb_state.clear();
      }      
  }    
 //  cout << "em_hmm_genotype::get_recomb_gradient() END" << endl;   
   return df_r_m;
}

void em_hmm_genotype::resolve_phase()
{
  cout << "em_hmm_genotype::resolve_phase() START" << endl;
  Phase PhaseObj;
  vector<vector<vector<int>>> orderedStates;
  vector<vector<vector<int>>> final_Phase;
  m_hap_data.m_theta.clear();
  m_hap_data.m_alpha.clear();
  m_current_states.clear();
  m_hap_data.m_recombinations.clear();
  
  //Read the statespace, model parameters for the instance with best parameters
  updateConvergedState();
  cout << "updateConvergedState() Done " << endl;
  //test_get_params();
   update_HMM_param();
   cout << "update_HMM_param() Done " << endl;
 //PhaseObj.resolvePhase(m_fwd_probs,m_clst_given_geno_v,m_hap_data,m_total_states,m_current_states,orderedStates,final_Phase);
  PhaseObj.resolvePhase(m_fwd_probs,m_clst_given_geno_v,m_hap_data,m_total_states,orderedStates,final_Phase);
  //m_params_file.close(); 
  cout << " PhaseObj.resolvePhase() Done " << endl;
  log_results(orderedStates,final_Phase);
 cout << "em_hmm_genotype::resolve_phase() END" << endl;
}

void em_hmm_genotype::log_results(vector<vector<vector<int>>> &orderedStates, vector<vector<vector<int>>> &final_Phase)
{
   cout << "em_hmm_genotype::log_results() START" << endl;
   ofstream clusters_file,phase_file; 
   int ind=0;
   clusters_file.open("Phase_clusters.txt");
   phase_file.open("Phased_haplo.txt");
   vector<vector<int>> indData;
   
  clusters_file << n_individuals<<endl;
  phase_file << n_individuals<<endl;
  clusters_file << n_markers<<endl;
  phase_file << n_markers << endl;
 
   for(int indCount=0;indCount<n_individuals;++indCount)
  {
      clusters_file << "# " << m_hap_data.m_input.m_ind_ID[indCount] << endl;
    
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
   cout << "em_hmm_genotype::log_results() END" << endl;
}

void em_hmm_genotype::updateConvergedState()
{
  // cout << "em_hmm_genotype::updateConvergedState() START" << endl;
   std::ifstream myReadFile;
   myReadFile.open("dakota_param.out");
   std::string line,lineObj,line1,line2;
   Double temp_Val,obj_function;
   istringstream buf_stream;
   int tempNum;
   std::vector<vector<int>> genoDataTemp;
   
   
   //################################## Read output file from Dakota and go to optimal params in that file ################################################3
   while(myReadFile.good())
   {
	getline(myReadFile,line); // get line from file
	   
	if(line.find("<<<<< Best parameters")!=string::npos) // search
	{
             break;
	}
    }
    
    for(int iCount=0;iCount < n_markers; iCount++)
    {
	 vector<Double> tempVec;
	 for(int jCount=0; jCount < n_clusters ;++jCount)
	 {
	      getline(myReadFile,line);	
	      line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
	      line = line.substr(0,line.find_first_of(" "));
	      temp_Val = stold(line);
	      tempVec.emplace_back(temp_Val);	     
	  }
	  m_hap_data.m_theta.emplace_back(tempVec);
     }
  
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	  vector<Double> tempVec;
	  
	  for(int jCount=0; jCount < n_clusters ;++jCount)
	  {
	       getline(myReadFile,line); 
	       line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
	       line = line.substr(0,line.find_first_of(" "));
	       temp_Val = stold(line);
	       tempVec.emplace_back(temp_Val);	          
	  }
	  m_hap_data.m_alpha.emplace_back(tempVec);
      }
   
      for(int iCount=0;iCount < n_markers-1; iCount++)
      {
	    getline(myReadFile,line);
	    line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
	    line = line.substr(0,line.find_first_of(" "));
	    temp_Val = stold(line);
	    m_hap_data.m_recombinations.emplace_back(temp_Val);	           
      }
      
   /*   getline(myReadFile,line);
      if(line.find("<<<<< Best objective function")!=string::npos) // search
      {
	  getline(myReadFile,line); 
	  line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::bind1st(std::not_equal_to<char>(), ' ')));
	  lineObj = line.substr(0,line.find_first_of(" "));
      }
    
     myReadFile.seekg (0, myReadFile.beg);
      
      while(myReadFile.good())	
      {
	getline(myReadFile,line); // get line from file
	if(line.find(lineObj)!=string::npos)// search
	{
	   getline(myReadFile,lineObj);
	   break;
	 }
	  line2 = line1;
	  line1 = line;
      }
      myReadFile.close();
     
      line2 = line2.substr(line2.find_last_of(" ")+1,(line2.length()));
      
      tempNum = stoi(line2);
      line2 = "statespace " + std::to_string(tempNum-1);
     
//#################################### Read the statespace for required instance from set of all recorded state spaces######################################
      myReadFile.open("state_space.txt");
      while(myReadFile.good())
      {
	getline(myReadFile,line); // get line from file
	
	if(line.find(line2)!=string::npos) // search
	{
	   break;
	}
      }
    
      m_current_states.clear();
      for(int icount=0;icount<n_individuals;++icount)
      {
	 getline(myReadFile,line);
	 getline(myReadFile,line);
	 vector<vector<int>> indData;
	
	 for(int jcount=0;jcount<n_markers;++jcount)
	 {
	   getline(myReadFile,line);
           buf_stream.str(line);
           genoDataTemp.emplace_back(istream_iterator<int>(buf_stream), istream_iterator<int>());
	   buf_stream.clear();   
	 }
	  
	 for(int j=0;j<n_markers;++j)
	 {
	    vector<int> temp;
	    for(int i =0; i < num_states;++i)
	      temp.emplace_back(genoDataTemp[j][i]);	      
	    indData.emplace_back(temp);
	 }
	  
	 genoDataTemp.clear();
	 m_current_states.push_back(indData);
      }*/
      myReadFile.close();
 
  //cout << "em_hmm_genotype::updateConvergedState() END" << endl;
}

void em_hmm_genotype::test_get_params()
{
      m_hap_data.compute_trans_prob_bet_clst();
      compute_jump_prob();  
    
      std::ifstream myReadFile; 
      string Strbuffer;
      istringstream buf_stream;
      vector<string> vData;
      
      myReadFile.open("CHKING.txt");  
      getline(myReadFile,Strbuffer);
      getline(myReadFile,Strbuffer);
      getline(myReadFile,Strbuffer);
      
      m_clst_given_v.clear();
      vector< double> dtemp_vec;
      for(int mrkCnt=0;mrkCnt<n_markers;++mrkCnt)
      {
	    getline(myReadFile,Strbuffer);
	    //cout << Strbuffer << endl;
	    getline(myReadFile,Strbuffer);
	    //cout << Strbuffer << endl;
	    buf_stream.str(Strbuffer);
	    dtemp_vec =  vector< double>(istream_iterator< double>(buf_stream), istream_iterator< double>()); 
	    vector<Double> markData;
	    for(auto &kv:dtemp_vec)
	    {
		  markData.emplace_back(kv);
	    }
	    m_clst_given_v.emplace_back(markData);
	    buf_stream.clear();
	}
      
     m_geno_given_clst_v.clear();
      getline(myReadFile,Strbuffer);
     
      for(int indCnt=0;indCnt<n_individuals;++indCnt)   
      {
	  getline(myReadFile,Strbuffer);
	  vector<vector<Double>> IndData;
	  for(int mrkCnt=0;mrkCnt<n_markers;++mrkCnt)
	  {
		 getline(myReadFile,Strbuffer);
		 getline(myReadFile,Strbuffer);
		 buf_stream.str(Strbuffer);
		 dtemp_vec =  vector< double>(istream_iterator< double>(buf_stream), istream_iterator< double>()); 
		 vector<Double> markData;
		  for(auto &kv:dtemp_vec)
		  {
			markData.emplace_back(kv);
		  }
		 IndData.emplace_back(markData);
		 buf_stream.clear();
	  }
	  m_geno_given_clst_v.emplace_back(IndData);
      }
      
      m_fwd_probs.clear();
      getline(myReadFile,Strbuffer);
      getline(myReadFile,Strbuffer);
     
      for(int indCnt=0;indCnt<n_individuals;++indCnt)   
      {
	  getline(myReadFile,Strbuffer);
	  vector<vector<Double>> IndData;
	  for(int mrkCnt=0;mrkCnt<n_markers;++mrkCnt)
	  {
		 getline(myReadFile,Strbuffer);
		 getline(myReadFile,Strbuffer);
		 buf_stream.str(Strbuffer);
		 dtemp_vec =  vector< double>(istream_iterator< double>(buf_stream), istream_iterator< double>()); 
		 vector<Double> markData;
		  for(auto &kv:dtemp_vec)
		  {
			markData.emplace_back(kv);
		  }
		 IndData.emplace_back(markData);
		 buf_stream.clear();
	  }
	  m_fwd_probs.emplace_back(IndData);
      }
      
       m_bckwd_probs.clear();
      getline(myReadFile,Strbuffer);
     getline(myReadFile,Strbuffer);
      
      for(int indCnt=0;indCnt<n_individuals;++indCnt)   
      {
	  getline(myReadFile,Strbuffer);
	  vector<vector<Double>> IndData;
	  for(int mrkCnt=0;mrkCnt<n_markers;++mrkCnt)
	  {
		 getline(myReadFile,Strbuffer);
		 getline(myReadFile,Strbuffer);
		 buf_stream.str(Strbuffer);
		 dtemp_vec =  vector< double>(istream_iterator< double>(buf_stream), istream_iterator< double>()); 
		 vector<Double> markData;
		  for(auto &kv:dtemp_vec)
		  {
			markData.emplace_back(kv);
		  }
		 buf_stream.clear();
	  }
	  m_bckwd_probs.emplace_back(IndData);
      }
      
      m_clst_given_geno_v.clear();
      getline(myReadFile,Strbuffer);
     getline(myReadFile,Strbuffer);
      
      for(int indCnt=0;indCnt<n_individuals;++indCnt)   
      {
	  getline(myReadFile,Strbuffer);
	  vector<vector<Double>> IndData;
	  for(int mrkCnt=0;mrkCnt<n_markers;++mrkCnt)
	  {
		 getline(myReadFile,Strbuffer);
		 getline(myReadFile,Strbuffer);
		 buf_stream.str(Strbuffer);
		  dtemp_vec =  vector< double>(istream_iterator< double>(buf_stream), istream_iterator< double>()); 
		 vector<Double> markData;
		  for(auto &kv:dtemp_vec)
		  {
			markData.emplace_back(kv);
		  }
		 IndData.emplace_back(markData);
		 buf_stream.clear();
	  }
	  m_clst_given_geno_v.emplace_back(IndData);
      }
      
   //   log_param();
}
