#include "globals.h"

Double dThreshold = pow(10,-30);
Double scaling_factor = pow(10,10);
double Recomb_InitialValue;
int n_individuals;
int n_markers;
int n_clusters;
int n_ploidy;
int num_states;
int n_scaler =0;
int n_runs=1;
int n_iterations =25;
int n_transitions = 12;
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

long int binomialCoef(int N,int r)
{
  long int iCoefficient=1,iResult;
  
  for(int iCount=0;iCount < r; iCount++)
  {
      iCoefficient *= (N-iCount) ;
     //iCoefficient = exp(log(iCoefficient)+log(N-iCount));
  }
  iResult = iCoefficient/factorial(r);
  //iResult = exp(log(iCoefficient)-log(factorial(r)));
  return(iResult) ;
  
}

Double binomialProb(int N,int r,Double p)
{
  long int iCoefficient = binomialCoef( N, r);
  Double dResult;
  dResult = iCoefficient*pow(p,r)*pow((1-p),(N-r));
  //dResult = exp(log(iCoefficient)+log(pow(p,r))+log(pow((1-p),(N-r))));
  return(dResult) ;
}

long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector< vector<int> >  &T)
{
  int i;
  long count = 0;
  
  if (n_chosen == len) 
  {
    vector<int> single;
    
    for (i = 0; i < len; i++)
    {
        single.emplace_back(input[got[i]]) ;
    }
    
    T.push_back(single);
    
    return 1;
  }
 
  for (i = at; i < max_types; i++) 
  {
    got[n_chosen] = i;
    count += choose(input,got, n_chosen + 1, len, i, max_types,T);
  }
  return count;
}

void vector_combination(vector<int >&input,vector<int> combination, int offset, int k,vector< vector<int>> &T) 
{
  if (k == 0) 
  {
    vector<int> single(combination);
    if(find( T.begin(),  T.end(), single) == T.end())
    T.push_back(single);
    return;
  }
  
  for (int i = offset; i <=  input.size() - k; ++i) 
  {
    combination.push_back(input[i]);
    vector_combination(input,combination,i+1, k-1,T);
    combination.pop_back();
  }
  
}

void vector_permutation(std::vector<int>& now,
std::vector<int> next,std::vector<vector<int>> &permutations)
{
  int size=now.size();
  
  if(size>0)
  {
	for(int cnt=0; cnt<size;++cnt)
	{
	      std::vector<int> vt;
	      std::vector<int>::const_iterator it=now.begin();
	      
	      for(int cnt1=0;cnt1<size;++cnt1)
	      {
		    if(cnt1==cnt)
		    {
		      ++it;
		      continue;
		    }
		    else
		   {
		      vt.push_back(*it);
		    }
		      ++it;
	      }
	      
	      std::vector<int>::const_iterator it1=now.begin();
	      --it1;
	      for(int cnt2=0;cnt2<=cnt;++cnt2)
	      {
		++it1;
	      }
	      next.push_back(*it1);
	      vector_permutation(vt,next,permutations);//,func);
	      next.pop_back();
        }
  }
  else
  {
	vector<int> single(next.begin(),next.end());
	if(std::find(permutations.begin(), permutations.end(), single)==permutations.end())
	{
	      permutations.push_back(single);
	}
  }
}
void gen_combi_states(vector< vector<int> > &T)
{
 
    vector<int> input(n_clusters); 
      
      for(int i= 0; i<n_clusters;i++)
      {
	input[i]= i;
      }
      
      vector<int> data(n_ploidy);
      choose(input,data, 0, n_ploidy, 0, n_clusters,T);  
    
      //gen_weights_states(T);
      return;
}
/*void gen_combi_states(vector< states > &T)
{
 
 vector<int> input(n_clusters); 
  
  for(int i= 0; i<n_clusters;i++)
  {
    input[i]= i;
  }
  
  vector<int> data(n_ploidy);
  choose(input,data, 0, n_ploidy, 0, n_clusters,T);  
 // gen_weights_states(nclusters, nploidy, T ,indiprobs);  
 gen_weights_states(T);
  return;
}*/

/*void sample(vector<states> &input,  vector<int> &exclude,vector<int> &output,int size_output)
//{
 
 //initialize random seed: 

  //srand (time(NULL));
  
  srand(time(0)+clock()+random());
  
  int j,k,choose,outputSize;
  int range = input.size();
  bool bErasebegin = false;
  outputSize = size_output;
 
  if(output.size()== 0)
  {
    output.push_back(-1); 
    outputSize = outputSize+1;
    bErasebegin = true;
  }
  
  while (output.size() < outputSize)
  {
    j = rand() % (range);       
    k = rand() % (range);
    choose = (input[j].weight < input[k].weight) ? k : j; 
    
    if( std::find(output.begin(), output.end(), choose)==output.end())
    {
      if( std::find(exclude.begin(), exclude.end(), choose)==exclude.end())
      {
	output.push_back(choose);    
      }
    }
    
  } 
  
  if(bErasebegin)
  {
    output.erase(output.begin()); 
  }
  
}*/

/*void sample(vector<states> &input,vector<int> &output,int size_output)
{
  vector<int> exclude{-1};
  sample(input,exclude,output,size_output);
}*/

/*void gen_weights_states( vector<states> &T)
{
  
      long int numerator_1;
      int loop_end= T.size();
      int k_count=0,value_count;
      Double value_temp;
      Double indprobs = (1.0/n_clusters) ;  
      numerator_1 = factorial(n_ploidy);
      
      for(int loop_count=0;  loop_count< loop_end; ++loop_count)
      {
	    k_count = 0;
	    value_temp = 1.0;
	    while(k_count < n_clusters)
	    {
		value_count = std::count(T[loop_count].values.begin(),T[loop_count].values.end(),k_count);
		value_temp = value_temp * pow(indprobs,value_count)/factorial(value_count);
		++k_count;
	    }
	    
	    value_temp *= numerator_1;
	    T[loop_count].weight = value_temp;	  
      }
}*/

int compare(const vector<int>& left, const vector<int>& right,vector<int>& out) 
{
      auto leftIt = left.begin();
      auto rightIt = right.begin();
      auto diff = 0;
      int pos=0;
      
      while (leftIt != left.end() && rightIt != right.end())
      {
	    if (*leftIt != *rightIt)
	    {
	      out.push_back(*rightIt);
	      diff++;
	    }
	    leftIt++;
	    rightIt++;
	    pos++;
      }

    return diff;
}

bool small_vectors(const vector<int> &a,const vector<int> &b) 
{
   return a.size() < b.size();
}

void get_map_vectors(vector<vector<int>>&input1, vector<int> &input2, vector<vector<std::pair<int,int>>> &output)
{
   // map the input1 and input1 in a unique way,  for eg input1[0]={0, 0, 1, 1} state={2,3,3,4} map= {(0, 2),(0, 3),(1, 3),(1, 4)}
      for(vector<vector<int>>::iterator haploCount= input1.begin();haploCount !=input1.end();++haploCount)
      {
	    vector<std::pair<int,int>> my_map; // new sub container      
	    std::transform((*haploCount).begin(), (*haploCount).end(), input2.begin(), 
	    std::inserter(my_map, my_map.end()), std::make_pair<int const&,int const&>);// push in a map for a single haplo vector and the state vector
			    
	    //std::sort(my_map.begin(),my_map.end(),pairCompare); // sort the map      
	    output.emplace_back(my_map);
	}
    
    //eliminate the duplicates
      std::sort(output.begin(),output.end()); // sort the map    
      output.erase(unique(output.begin(),output.end()), output.end()); 
}

Double  randZeroToOne()
{
    srand((unsigned)time(NULL));
    return rand() / (RAND_MAX + 1.);
}

/*void print_values( int size,  const gsl_vector *x)
{
   cout << "print_values " << endl;
    for(int i=0;i<size;++i)
    {
	cout << gsl_vector_get (x,i) << " ";
    }
   cout << endl;
}*/