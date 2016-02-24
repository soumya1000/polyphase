#include "parsedata.h"

 // All global variables
int n_individuals;
int n_markers;
int n_clusters;
int n_ploidy;
Double Recomb_InitialValue;
int num_states;

data_genotypes::data_genotypes()
{
  //cout << "data_genotypes::data_genotypes()" << endl;
}

data_genotypes::~data_genotypes()
{
}

void data_genotypes::read_paramsFile()
{
      string line;
      std::ifstream paramsFile ("Poly_phase_params.txt");
      
      try
     {
	    if (paramsFile.is_open()) 
	    {
		  //getline(paramsFile,line);
		   getline(paramsFile,line);
		  
		  getline(paramsFile,line);		  
		  line = line.substr(line.find_first_of(":"),line.length());
	          line = line.substr(line.find_first_of(" "),line.length());		
		  n_clusters = stoi(line);
		  
		  getline(paramsFile,line);
		  line = line.substr(line.find_first_of(":"),line.length());
	          line = line.substr(line.find_first_of(" "),line.length());			
		  Recomb_InitialValue = stold(line);
		  
		  getline(paramsFile,line);		  
		  line = line.substr(line.find_first_of(":"),line.length());
	          line = line.substr(line.find_first_of(" "),line.length());
		 
		  num_states = stoi(line);
	    }
	    else
	    {
		  cout << "Please provide params in Poly_phase_params.txt" << endl;
		  exit(1);
	    }
     }
     catch(const std::exception& e)
    {
         cout << "Error in reading Poly_phase_params.txt " ;
	std::cout << e.what() << endl;
	exit(1);
    }
}

void data_genotypes::parse_input(string iFilename)
{
    //cout << "data_genotypes::parse_input() START" << endl; 
      read_paramsFile();
      std::ifstream myReadFile (iFilename); 
      string Strbuffer;
      istringstream buf_stream;
     
      std::vector<vector<int>> genoDataTemp;
      if (myReadFile.is_open()) 
      {
	    getline(myReadFile,Strbuffer);
	    n_ploidy = std::stoi(Strbuffer) ;
	    
	    getline(myReadFile,Strbuffer);
	    n_individuals = std::stoi(Strbuffer) ;
	    
	    getline(myReadFile,Strbuffer);
	    n_markers = std::stoi(Strbuffer) ;
	    
	    getline(myReadFile,Strbuffer);
	    Strbuffer.erase(Strbuffer.begin(),Strbuffer.begin()+1);
	    
	    buf_stream.str(Strbuffer);
	    m_physical_positions = vector<int>(istream_iterator<int>(buf_stream), istream_iterator<int>()); 
	    for(int i=0;i<m_physical_positions.size()-1;++i)
	    {
		  m_physical_distances.push_back(m_physical_positions[i+1]-m_physical_positions[i]);
	    }
	    buf_stream.clear();
	    while (!myReadFile.eof()) 
	    {
		  getline(myReadFile,Strbuffer);
		  buf_stream.str(Strbuffer);
		  std::vector<vector<int>>genoData;
		  if(Strbuffer.compare(0,1,"#")==0)
		  {
			  m_ind_ID.emplace_back(Strbuffer.substr(Strbuffer.find_first_of(" "),Strbuffer.length()));
			  buf_stream.clear();			  
			  //read the unphased haplotypes and store them
			  for(int hapCount=0;hapCount<n_ploidy;++hapCount)
			  {
				getline(myReadFile,Strbuffer);
				buf_stream.str(Strbuffer);
				genoDataTemp.emplace_back(istream_iterator<int>(buf_stream), istream_iterator<int>());
				buf_stream.clear();
				/* std::transform(genoData.begin(), genoData.end(), genoDataTemp.begin(),
				  genoData.begin(), std::plus<int>());*/
			  }			
			  for(int j =0; j < n_markers;++j) 
			  {
				vector<int> temp;
				for(int i =0; i < genoDataTemp.size();++i)
				    temp.emplace_back(genoDataTemp[i][j]);
				  genoData.emplace_back(temp);
			  }
			  m_genotypes.push_back(genoData);
			  genoDataTemp.clear();
		    }
	      }
	      myReadFile.close();
	}
	else
	{
	    cout << "Please provide input data in the file,"<< iFilename << endl;
	    exit(1);
	}
      
      
      //cout << "data_genotypes::parse_input() END" << endl; 
}


void data_genotypes::log_variables()
{
      std::ofstream file_knownParams;
      file_knownParams.open("phase_param.in"); 
      file_knownParams << n_ploidy << endl <<  n_individuals << endl << n_markers << endl<<n_clusters <<endl << num_states << endl ;
      
      for(const auto& kvp : m_physical_positions) 
      {
	    file_knownParams << kvp << "\t";
      }
      
      int ind=0; 
      file_knownParams << endl;
      for(const auto& kvp : m_genotypes) 
      {
	    //cout << "ind" << ind << " ";     
	    for(auto kv : kvp)
	    {
		for(auto v : kv)
		  file_knownParams << v<<" ";	
		file_knownParams << endl;
	    }
      }
      
      for(const auto& kvp : m_ind_ID) 
      {
	    file_knownParams << kvp << "\t";
      }
      file_knownParams.close();
}

void data_genotypes::read_variables()
{
   std::ifstream myReadFile;
   std::string line;
   istringstream buf_stream;
   myReadFile.open("phase_param.in");
   istream_iterator<int> iter_obj;
   getline(myReadFile,line);	
   n_ploidy = std::stoi(line) ;
   
   getline(myReadFile,line);
   n_individuals = std::stoi(line) ;
    
   getline(myReadFile,line);
   n_markers = std::stoi(line) ;
   
   getline(myReadFile,line);
   n_clusters = std::stoi(line) ;
   
   getline(myReadFile,line);
   num_states = std::stoi(line) ;
   
   getline(myReadFile,line);
   buf_stream.str(line);
   m_physical_positions = vector<int>(istream_iterator<int>(buf_stream), istream_iterator<int>()); 
   
   for(int i=0;i<m_physical_positions.size()-1;++i)
   {
      m_physical_distances.push_back(m_physical_positions[i+1]-m_physical_positions[i]);
   }
   buf_stream.clear();
    
   for(int indCount=0;indCount<n_individuals;++indCount)
   {
	vector<vector<int>> indData;
	for(int markCount=0;markCount<n_markers;++markCount)
	{
	      vector<int> markData;
	      getline(myReadFile,line);
	      buf_stream.str(line);
	      markData=vector<int>(istream_iterator<int>(buf_stream), istream_iterator<int>()); 
	      buf_stream.clear();
	      indData.emplace_back(markData);
	}
	m_genotypes.emplace_back(indData);
   }
   getline(myReadFile,line);
   buf_stream.str(line);
   m_ind_ID= vector<string>(istream_iterator<string>(buf_stream), istream_iterator<string>());
   myReadFile.close();
   
  /* for(const auto& kvp : m_physical_positions) 
  {
     cout << kvp << "\t";
  }
  
  int ind=0; 
  cout << endl;
  for(const auto& kvp : m_genotypes) 
  {
    //cout << "ind" << ind << " ";     
    for(auto kv : kvp)
    {
	for(auto v : kv)
	  cout << v<<" ";	
	cout << endl;
    }
  }
  
  for(const auto& kvp : m_ind_ID) 
  {
     cout << kvp << "\t";
  }*/
}

//Test code
/*int main(int argc, char *argv[])
{
  data_genotypes inputObj;
  inputObj.parse_input();
  inputObj.print_variables();
}*/