#include "haplo.h"


//*************************************************************//
// CLASS Chaplotypes
//*************************************************************//
Chaplotypes::Chaplotypes()
{
  //cout << "Chaplotypes::Chaplotypes()" << endl;
 
  //m_params_file_haplos.open ("Model_param.txt"); 
}

Chaplotypes::~Chaplotypes()
{
  //m_params_file_haplos.close();
  //cout << "Chaplotypes::~Chaplotypes()" << endl;
  //print_param();
  //log_param_hap();
  //m_params_file_haplos.close();
}

void Chaplotypes::initialise_param()
{
	//initialise theta matrix
	std::random_device rd;
	std::mt19937 gen(rd());
	//std::default_random_engine generator;
	
	std::uniform_real_distribution<double> distribution(0.0001,0.9999);
	
	// THETA values
	m_theta.clear();
	for(int iCount=0;iCount < n_markers; ++iCount)
	{
	    vector<double> thetaTemp;
	  
	    for(int jCount=0; jCount < n_clusters ;++jCount)
	    {
	      thetaTemp.emplace_back(distribution(rd));        
	    }
	  
	    m_theta.emplace_back(thetaTemp);
	}
	
	//initialse cluster frequencies (ALPHA)  
	//LATER : assign initial values from dirichlet distribution
	 double normalise; 
	 m_alpha.clear();
        for(int iCount=0;iCount < n_markers; ++iCount)
	{  
	    vector<double> alphaTemp;  
	    normalise = 0.0;
	    for(int jCount=0;jCount< n_clusters;++jCount)
	    {
		alphaTemp.emplace_back(1.0/n_clusters); 
		//alphaTemp.emplace_back(distribution(rd)); 
		normalise = normalise + alphaTemp[jCount];
	     } 
	     for(int jCount=0; jCount < n_clusters ;++jCount)
	     {
	       alphaTemp[jCount] = alphaTemp[jCount]/normalise ;
	     }
		m_alpha.emplace_back(alphaTemp);
	  }   
          double dtemp;
	  m_recombinations.clear();
	  for(int iCount =0;iCount < (n_markers-1) ;++iCount)
	  {
	      dtemp= 1-exp(-1*Recomb_InitialValue*m_input.m_physical_distances[iCount]);
	      m_recombinations.emplace_back(dtemp);
	   }
	//print_param();
}

void Chaplotypes::compute_trans_prob_bet_clst()
{
   //cout << "Chaplotypes::compute_trans_prob_bet_clst() START" << endl;
    Double dTempJumpprob=0,dTemp;

      vector<double> alphaCurMarker;
      m_trans_prob_bet_clst.clear();
      int loopEnd = n_clusters*n_clusters;
      map<pair<int,int>,Double> tempData;
      vector <pair<int,int>> tempDataVec;
      for(int kCount=0; kCount < n_clusters; ++kCount)
      { 
	   for(int lCount=0; lCount < n_clusters; ++lCount)
	   {    
		tempData[std::make_pair (kCount,lCount)] =0.0 ;
		tempDataVec.emplace_back(std::make_pair (kCount,lCount));
	   }  
      } 

        for(int jCount=0; jCount < n_markers-1; jCount++)
	{ 
	    map<pair<int,int>,Double> markerData = tempData;
	    alphaCurMarker = m_alpha[jCount];
	    parallel_for(int(0),loopEnd,[this,&markerData,&jCount,&tempData,&alphaCurMarker,&tempDataVec](int iCounter)throw()
	    {
	         Double dTemp =  m_recombinations[jCount];
		 Double dTempJumpprob = (dTemp)*alphaCurMarker[tempDataVec[iCounter].second];   
		 
		 if(tempDataVec[iCounter].first == tempDataVec[iCounter].second)
		 {
		     dTempJumpprob = (dTempJumpprob + 1-dTemp);
		 }
	          markerData[ tempDataVec[iCounter]] = dTempJumpprob;
	    });
	  
	    m_trans_prob_bet_clst.emplace_back(markerData);      
	} 
    /* int markCount =0;
      for(auto &it:m_trans_prob_bet_clst) 
      {
	    cout << "markCount "<< markCount << endl;
	    for(auto &jt:it) 
	    {
		cout <<  jt.first.first << "," << jt.first.second << " : " << jt.second << endl;
	    }
	    ++markCount;
      }*/
  // cout << "Chaplotypes::compute_trans_prob_bet_clst() END" << endl;
  
}

Double Chaplotypes::compute_trans_prob_bet_clst_tuples(int frmMarker,vector<vector<int>>&fromstatePerms,vector<int>&toState)
{
   //cout << "Chaplotypes::compute_trans_prob_bet_clst_tuples START" << endl;
  
   map<pair<int,int>,Double> markerData;  
   Double transValue = 0.0,tempTrans;
    
   markerData = m_trans_prob_bet_clst[frmMarker];

   vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
   get_map_vectors(fromstatePerms, toState, Map_Valid_stateTuples);
   
   for(auto &kvp:Map_Valid_stateTuples)
   {
      tempTrans = 1.0;
      for(auto kv:kvp)
      {
	  //tempTrans = exp(log(tempTrans) + log(markerData[kv]));
	  tempTrans *= markerData[kv];
      }
      transValue +=  tempTrans;
   }
   
    return transValue;  
}

void  Chaplotypes::print_param(void)
{
   cout << "Theta  " << endl;    
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      cout << "Marker: "<< iCount << endl;
      for(int jCount=0; jCount <n_clusters ;jCount++)
      {
	  cout << m_theta[iCount][jCount]<< " ";
      }
      cout << endl;
   }
    
   cout << "alpha " << endl;   
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      cout << "Marker: "<< iCount << endl;
      for(int jCount=0; jCount <n_clusters ;jCount++)
      {
        cout << m_alpha[iCount][jCount] << " ";
      }
      cout << endl;
    }
    
    cout << "Recombination values :" << endl;
    for(int iCount =0;iCount < n_markers-1 ;iCount++)
    {
      cout << m_recombinations[iCount] << endl;
    }
}

void Chaplotypes::log_param_hap(string fName,bool recombDirect)
{
      ofstream params_file(fName,ios::out | ios::app);
      vector<double> recombs_final;
      //params_file << endl << "**********" << endl;
      
     params_file << "alpha " << endl;        
      for(auto &kvp:m_alpha)
      {
	  //params_file << "Marker: "<< ++marker << endl;
	  for(auto kv:kvp)
	  {
	    params_file  << kv<< " ";
	  }
	  params_file << endl;
	}
	
      params_file << "Theta  " << endl;  
      int marker =-1;
      double temp;
      for(auto &kvp:m_theta)
      {
	  //params_file << "Marker: "<< ++marker << endl;
	  for(auto kv:kvp)
	  {
	      params_file << kv<< " ";
	  }
	  params_file << endl;
      }
      marker =-1; 
      marker =0; 	
	for(auto &kvp:m_recombinations)
	{
	    temp = log(1/(1-kvp))/m_input.m_physical_distances[marker];
	    recombs_final.emplace_back(temp);
	    ++marker;
	}
	if(recombDirect==true)
	{
	    params_file << endl<<"Recombination values Final:" << endl;
	    for(auto &kvp:recombs_final)
	    {
		params_file << kvp << " ";
	    }
	}
	else
	{
	    params_file << "1-e^(-rmdm) :" << endl;
	    for(auto &kvp:m_recombinations)
	    {
		params_file << kvp << " ";
	    }
	  
	}
	params_file << endl << "**********" << endl;
	params_file.close();
}

void Chaplotypes::initialize(string ipFileName)
{
    m_input.parse_input(ipFileName);
    initialise_param();
}
