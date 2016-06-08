#include "haplo.h"


//*************************************************************//
// CLASS Chaplotypes
//*************************************************************//
Chaplotypes::Chaplotypes()
{
  cout << "Chaplotypes::Chaplotypes()" << endl;
 
  m_params_file_haplos.open ("Model_param.txt"); 
}

Chaplotypes::~Chaplotypes()
{
  //m_params_file_haplos.close();
  cout << "Chaplotypes::~Chaplotypes()" << endl;
  //print_param();
  //log_param_hap();
  m_params_file_haplos.close();
}

void Chaplotypes::initialise_param(string ipFileName)
{
	m_input.parse_input(ipFileName);

	//initialise theta matrix
	std::random_device rd;
	std::mt19937 gen(rd());
	//std::default_random_engine generator;
	
	std::uniform_real_distribution<double> distribution(0.0001,0.9999);
	
	// THETA values
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
        for(int iCount=0;iCount < n_markers; ++iCount)
	{  
	    vector<double> alphaTemp;  
	    normalise = 0.0;
	    for(int jCount=0;jCount< n_clusters;++jCount)
	    {
		//alphaTemp.emplace_back(1.0/n_clusters); 
		alphaTemp.emplace_back(distribution(rd)); 
		normalise = normalise + alphaTemp[jCount];
	     } 
	     for(int jCount=0; jCount < n_clusters ;++jCount)
	     {
	       alphaTemp[jCount] = alphaTemp[jCount]/normalise ;
	      }
		m_alpha.emplace_back(alphaTemp);
	  }   

	  for(int iCount =0;iCount < (n_markers-1) ;++iCount)
	  {
	      m_recombinations.emplace_back(Recomb_InitialValue);
	   }
	//print_param();
}

void Chaplotypes::compute_trans_prob_bet_clst()
{
   //cout << "Chaplotypes::compute_trans_prob_bet_clst() START" << endl;
    Double dTempJumpprob=0,dTemp;

      vector<double> alphaCurMarker;
      
	for(int jCount=0; jCount < n_markers-1; jCount++)
	{ 
	    map<pair<int,int>,Double> markerData;
	    alphaCurMarker = m_alpha[jCount];
	  
	    for(int kCount=0; kCount < n_clusters; kCount++)
	    { 
		for(int lCount=0; lCount < n_clusters; lCount++)
		{    
		    dTemp = exp(-1*m_recombinations[jCount]*m_input.m_physical_distances[jCount]);
		    //dTempJumpprob = (1-dTemp)*alphaCurMarker[lCount];    
		    dTempJumpprob = exp(log(1-dTemp)+ log(alphaCurMarker[lCount]));  
		    
		    if( kCount == lCount)
		    {
		      dTempJumpprob = dTempJumpprob + dTemp;
		      //dTempJumpprob = exp(add_log(log(dTempJumpprob),log(dTemp)));
		    }
	      
		    markerData[ std::make_pair (kCount,lCount)] = dTempJumpprob;
		}  
	    } 
	    m_trans_prob_bet_clst.emplace_back(markerData);      
	} 

      /*std::map<pair<int,int>, double> x;
      for(size_t i= 0;i < m_trans_prob_bet_clst.size();i++) 
      {
	x = m_trans_prob_bet_clst[i];
	for(  std::map<pair<int,int>, double>::iterator j= x.begin();j!=x.end();j++) 
	{
	  cout <<  (*j).first.first << "," << (*j).first.second << " : " << (*j).second << endl;
	}
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
	  tempTrans = exp(log(tempTrans) + log(markerData[kv]));
      }
      transValue = transValue + tempTrans;
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
	  cout << " " << m_theta[iCount][jCount];
      }
      cout << endl;
   }
    
   cout << "alpha " << endl;   
   for(int iCount=0;iCount < n_markers; iCount++)
   {
      cout << "Marker: "<< iCount << endl;
      for(int jCount=0; jCount <n_clusters ;jCount++)
      {
        cout << " " << m_alpha[iCount][jCount];
      }
      cout << endl;
    }
    
    cout << "Recombination values :" << endl;
    for(int iCount =0;iCount < n_markers-1 ;iCount++)
    {
      cout << m_recombinations[iCount] << endl;
    }
}

void Chaplotypes::log_param_hap(void)
{
   m_params_file_haplos << endl << "**********************************************************" << endl;
   m_params_file_haplos << "Theta  " << endl;  
   int marker =-1;
   for(auto &kvp:m_theta)
   {
      m_params_file_haplos << "Marker: "<< ++marker << endl;
      for(auto kv:kvp)
      {
	  m_params_file_haplos << " " << kv;
      }
      m_params_file_haplos << endl;
   }
   marker =-1; 
   m_params_file_haplos << "alpha " << endl;   
   
   for(auto &kvp:m_alpha)
   {
      m_params_file_haplos << "Marker: "<< ++marker << endl;
      for(auto kv:kvp)
      {
        m_params_file_haplos << " " << kv;
      }
      m_params_file_haplos << endl;
    }
    marker =-1; 
    m_params_file_haplos << "Recombination values :" << endl;
    for(auto &kvp:m_recombinations)
    {
      m_params_file_haplos << kvp << " ";
    }
     m_params_file_haplos << endl << "*******************************************************" << endl;
}
