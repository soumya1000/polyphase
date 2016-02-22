#include "haplo.h"

//*************************************************************//
// CLASS HAPLOTYPES
//*************************************************************//
Chaplotypes::Chaplotypes()
{
  
  //cout << "Chaplotypes::Chaplotypes()" << endl;
  
}

Chaplotypes::~Chaplotypes()
{
  ////cout << "Chaplotypes::~Chaplotypes()" << endl;
  //print_param();
  log_param_hap();
}

void Chaplotypes::initialise_param()
{
 // cout << "Chaplotypes::initialise_param() START" << endl;
  //initialise theta matrix
  std::random_device rd;
  std::mt19937 gen(rd());
  //std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0001,0.9999);
  m_theta.clear();
  m_alpha.clear();
  for(int iCount=0;iCount < n_markers; iCount++)
  {
      vector<Double> thetaTemp;
     
      for(int jCount=0; jCount < n_clusters ;jCount++)
      {
        thetaTemp.push_back(distribution(rd));        
      }
     
      m_theta.push_back(thetaTemp);
   }
      
   //initialse cluster frequencies (alpha) 
 
  //LATER : assign initial values from dirichlet distribution
 Double normalise;
 
  for(int iCount=0;iCount < n_markers; iCount++)
  {  
    vector<Double> alphaTemp;  
    normalise = 0.0;
    for(int jCount=0;jCount< n_clusters;jCount++)
    {
      //alphaTemp.push_back(1.0/n_clusters); 
      alphaTemp.emplace_back(distribution(rd)); 
      normalise = normalise + alphaTemp[jCount];
    } 
    for(int jCount=0; jCount < n_clusters ;jCount++)
    {
        alphaTemp[jCount] = alphaTemp[jCount]/normalise ;
    }
    m_alpha.emplace_back(alphaTemp);
  } 
  for(int iCount =0;iCount < n_markers-1 ;iCount++)
  {
    m_recombinations.push_back(Recomb_InitialValue);
  }
  
  //cout << "Chaplotypes::initialise_param() END" << endl;
}

void Chaplotypes::initial_call_actions(int num_clusters)
{
  //cout << "Chaplotypes::initial_call_actions() START" << endl;
   m_input.parse_input(); 
   //m_input.print_variables();
  //initialise clusters
  
  n_clusters = num_clusters; 
  m_input.log_variables();
  initialise_param();
  optimise_param_helper();
  std::remove("state_space.txt"); 
  std::remove("Hmm_params.txt");
  //cout << "Chaplotypes::initial_call_actions() END" << endl;
}
void Chaplotypes::restore_prevstate()
{
 
   m_input.read_variables();     
}
void Chaplotypes::compute_trans_prob_bet_clst()
{
  //cout << " CHK:  Chaplotypes::compute_trans_prob_bet_clst START" << endl;
  
 Double dTempJumpprob=0;
 Double dTemp,dTempval;
 
  vector<Double> alphaCurMarker;
  
   for(int jCount=0; jCount < n_markers-1; jCount++)
    { 
      map<pair<int,int>,Double> markerData ;
 
      alphaCurMarker = m_alpha[jCount];
       
	for(int kCount=0; kCount < n_clusters; kCount++)
	{ 
	    for(int lCount=0; lCount < n_clusters; lCount++)
	    {    
	      
		dTemp = exp(-1.0*m_input.m_physical_distances[jCount] *m_recombinations[jCount]);
		dTempJumpprob = (1-dTemp)*alphaCurMarker[lCount];    
	      		
		if( kCount == lCount)
		{
		  dTempJumpprob = dTempJumpprob + dTemp;
		}
           
		markerData[ std::make_pair (kCount,lCount)] = dTempJumpprob;
	    }  
	} 
	m_trans_prob_bet_clst.push_back(markerData);    
	
    } 

   // cout << " CHK:  Chaplotypes::compute_trans_prob_bet_clst END" << endl;
  /*std::map<pair<int,int>, long double> x;
  for(size_t i= 0;i < m_trans_prob_bet_clst.size();i++) 
  {
    x = m_trans_prob_bet_clst[i];
    for(  std::map<pair<int,int>, long double>::iterator j= x.begin();j!=x.end();j++) 
    {
      cout <<  (*j).first.first << "," << (*j).first.second << " : " << (*j).second << endl;
    }
  }*/
}

Double Chaplotypes::compute_trans_prob_bet_clst_tuples(int frmMarker,vector<vector<int>>&fromstatePerms,vector<int>&toState)
{
   //cout << "Chaplotypes::compute_trans_prob_bet_clst_tuples START" << endl;
   
   map<pair<int,int>,Double> markerData;  
   Double tempTrans,transValue=0,print=0;
   vector<int> frmstatevec;
   markerData = m_trans_prob_bet_clst[frmMarker];

   vector< vector<std::pair<int,int>>> Map_Valid_stateTuples;
   get_map_vectors(fromstatePerms, toState, Map_Valid_stateTuples);
   
  for(auto &kvp:Map_Valid_stateTuples)
   {
      tempTrans = 1.0;
      for(auto kv:kvp)
      {
	 //mpfr_mul(tempTrans.backend().data(), tempTrans.backend().data(),  markerData[kv].backend().data(), GMP_RNDD);
	  tempTrans *= markerData[kv];	 
      }
      //
       transValue = transValue + tempTrans;
      //mpfr_add (transValue.backend().data(), transValue.backend().data(), tempTrans.backend().data(), GMP_RNDU);
     
   }

  //  cout << "Chaplotypes::compute_trans_prob_bet_clst_tuples END" << endl;    
   return  transValue; 
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
   std::string line;
   std::ifstream myReadFile;
   myReadFile.open("mod_param_Haplo.txt"); 
   
   int count =0;
   while(myReadFile)
   {
      getline(myReadFile,line);
      
      if(line.find("Theta")!=string::npos) // search
	  ++count;
   } 
   myReadFile.close(); 
  
  /* m_params_file_haplos << endl << "**********************************************************" << endl;
   m_params_file_haplos << "Theta  " << count <<endl;  
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
     m_params_file_haplos << endl << "*******************************************************" << endl;*/
}

void Chaplotypes::optimise_param_helper()
{
  // cout<< "CHK em_hmm_genotype::optimise_param_helper() START" << endl;
  //ofstream input_file; 
   std::ofstream input_file("dakota_param.in");
 
   int ind=0;
   int  n_max_iterations = 15;//25
   Double d_convergence_tolerance = 1e-4;//1e-10; 
   Double  step_size = 1e-3;//1e-2;
   Double upper_bound = 0.99999;
   Double lower_bound = 1e-10;
   Double upper_bound_recomb = 0.1;
    Double lower_bound_recomb = 1e-10;
   int num_Params;  
   num_Params =  (n_markers* n_clusters * 2) + (n_markers-1) ;//theta, alpha, r .  
   string var_name;
   string method_name = "conmin_frcg";
   
   
      //input_file.open("phase_param.in");
  
      input_file << "environment" << endl;
      input_file << " tabular_graphics_data" << endl;
      input_file << "  tabular_graphics_file = 'poly_phase.dat' " << endl;
      input_file << endl <<"method" <<endl;
      
      input_file << " max_iterations = " << n_max_iterations << endl;
      
      input_file << " convergence_tolerance =  " << d_convergence_tolerance << endl;
      input_file << method_name << endl;
      input_file << endl <<"model" << endl;
      input_file << " single" << endl;
      input_file << endl <<"variables" << endl;
      input_file <<  " continuous_design = " << num_Params << endl;   
      input_file << "initial_point " ;//   -1.2      1.0
    
     
      //list the initial values of the parameters
    
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	  for(int jCount=0; jCount < n_clusters ;jCount++)
	  {
	    input_file << "\t" << m_theta[iCount][jCount];
	  }
       }
    
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	  for(int jCount=0; jCount < n_clusters; jCount++)
	  {
	      input_file << "\t" << m_alpha[iCount][jCount];

	  }
      }
     for(int iCount=0;iCount < n_markers-1; iCount++)
       {
	    input_file << "\t" << m_recombinations[iCount];

       }
   
      input_file << endl << "lower_bounds" ;
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	  for(int jCount=0; jCount < n_clusters ;jCount++)
	  {
	    input_file << "\t" << lower_bound << "\t" << lower_bound;
	  }
	
      }
      
      for(int iCount=0;iCount < n_markers-1; iCount++)
      {
	  input_file << "\t" << lower_bound_recomb;
      }
     
      input_file << endl << "upper_bounds" ; 
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	  for(int jCount=0; jCount < n_clusters ;jCount++)
	  {
	      input_file << "\t" << upper_bound << "\t" << upper_bound;
	  }
      }
      for(int iCount=0;iCount < n_markers-1; iCount++)
      {
        input_file << "\t" << upper_bound_recomb;
      }
      input_file << endl << "descriptors" ; 
      
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	for(int jCount=0; jCount < n_clusters ;jCount++)
	{
	  var_name =  std::to_string(iCount)+ std::to_string(jCount);
	  input_file << "\t" << "'theta_" + var_name + "'" ;
	}
      }
      for(int iCount=0;iCount < n_markers; iCount++)
      {
	for(int jCount=0; jCount < n_clusters ;jCount++)
	{
	   var_name =  std::to_string(iCount)+ std::to_string(jCount);
	   input_file << "\t" << "'alpha_"+ var_name + "'";
	}
      }
      for(int iCount=0;iCount < n_markers-1; iCount++)
      {
	  input_file << "\t" << "'recomb_"+ std::to_string(iCount)+"'" ;
      }
      
      
    /* input_file << endl << endl << "interface" << endl << " analysis_driver = 'phase_sim'" << endl << " file_tag" << endl 
    << " file_save" <<endl << " fork" << endl<< " parameters_file = 'params.in'" << endl << "results_file    = 'results.out'" << endl;*/
        
     input_file << endl << endl << "interface" << endl << " analysis_driver = 'phase_sim'" << endl 
     << " fork" << endl<< " parameters_file = 'params.in'" << endl << "results_file    = 'results.out'" << endl;
      
      input_file << endl << "responses" << endl << " objective_functions = 1" << endl;
    
      input_file << "numerical_gradients" << endl; 

      input_file <<"method_source dakota" << endl<< " interval_type forward" << endl << " fd_gradient_step_size = " << step_size
       << endl << " no_hessians" ;
   
  input_file.close();
 //cout<< "CHK em_hmm_genotype::optimise_param_helper() END" << endl;  
   
} 
  