#include "globals.h"


  Double dThreshold = pow(10,-10);
  Double scaling_factor = pow(10,10);
 


int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

long int binomialCoef(int N,int r)
{
  long int iCoefficient=1,iResult;
  
  for(int iCount=0;iCount < r; iCount++)
  {
    //iCoefficient = iCoefficient*(N-iCount) ;
     iCoefficient = exp(log(iCoefficient)+log(N-iCount));
  }
  //iResult = iCoefficient/factorial(r);
  iResult = exp(log(iCoefficient)-log(factorial(r)));
  return(iResult) ;  
}

Double binomialProb(int N,int r,Double p)
{
  long int iCoefficient = binomialCoef( N, r);
  Double dResult;
  //dResult = iCoefficient*pow(p,r)*pow((1-p),(N-r));
  dResult = exp(log(iCoefficient)+log(pow(p,r))+log(pow((1-p),(N-r))));
  return(dResult) ;
}

long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector< states > &T)
{
  int i;
  long count = 0;
  
  if (n_chosen == len) 
  {
    states single;
    
    for (i = 0; i < len; i++)
    {
        single.values.push_back(input[got[i]]) ;
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
long choose(vector<int > &input,vector<int > &got, int n_chosen, int len, int at, int max_types,vector< vector<int>> &T)
{
  int i;
  long count = 0;
  
  if (n_chosen == len) 
  {
    vector<int> single;
    
    for (i = 0; i < len; i++)
    {
        single.push_back(input[got[i]]) ;
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
void gen_combi_states(int nploidy,int nclusters,vector< vector<int> > &T)
{
      vector<int> input(nclusters); 
  
      for(int i= 0; i<nclusters;i++)
      {
	input[i]= i;
      }
      vector<int> data(nploidy);
      choose(input,data, 0, nploidy, 0, nclusters,T);
}
/*void gen_combi_states(int nploidy,int nclusters,vector< states > &T)
{
 
 vector<int> input(nclusters); 
  
  for(int i= 0; i<nclusters;i++)
  {
    input[i]= i;
  }
  
  vector<int> data(nploidy);
  choose(input,data, 0, nploidy, 0, nclusters,T);
  
  vector<Double> indiprobs;
  for(int i= 0; i<nclusters;i++)
  {
    indiprobs.push_back( 1.0/nclusters);
  }
  
  gen_weights_states(nclusters, nploidy, T ,indiprobs);  
  return;
}*/


void sample(vector<states> &input,  vector<int> &exclude,vector<int> &output,int size_output)
{
 
 /*initialize random seed: */

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
  
}
void sample(vector<states> &input,vector<int> &output,int size_output)
{
  vector<int> exclude{-1};
  sample(input,exclude,output,size_output);
}
void gen_weights_states(int n, int k, vector<states> &T, vector<Double> &indiprobs)
{
  
  int numerator_1;
  int loop_count=0,k_count=0,value_count;
  Double value_temp;
  vector<int> k_vals;
  
  while(k_count < k)
  {
     k_vals.push_back(loop_count);
     k_count++;
  }
  
  loop_count= 0;   
  numerator_1 = factorial(n);
  
  while(loop_count < T.size())
  {
    k_count = 0;
    value_temp = 1.0;
    
    while(k_count < k)
    {
      value_count = std::count(T[loop_count].values.begin(),T[loop_count].values.end(),k_count);
      value_temp = value_temp * pow(indiprobs[k_count],value_count)/factorial(value_count);
       //value_temp = exp(log(value_temp)+log(pow(indiprobs[k_count],value_count))-log(factorial(value_count)));
      k_count++;
    }
    
    //value_temp = pow(10,(log10(value_temp)+log10(numerator_1)));
    value_temp = value_temp*numerator_1;
    T[loop_count].weight = value_temp;
    loop_count++;
  }
}

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


void normalize_vector(vector<double> &inout)
{
  double sum_elements = std::accumulate(inout.begin(), inout.end(),0.0);
  
  if(sum_elements <= 1.0)
    return;
 
   for(auto  &it:inout)
     it = it/sum_elements;
}

// projects the ind values of vector between 0 and 1
void project_range0_1(vector<double> &inout)
{
  long int num_copy, num_ClosestVal;
  
  for(auto &it:inout)
  {
      if(abs(it)>1)
      {
	 it = atan(it) ;
      }    
  }  
}
void get_map_vectors(vector<vector<int>>&input1, vector<int> &input2, vector<vector<std::pair<int,int>>> &output)
{
   // map the input1 and input1 in a unique way,  for eg input1[0]={0, 0, 1, 1} state={2,3,3,4} map= {(0, 2),(0, 3),(1, 3),(1, 4)}
  for(vector<vector<int>>::iterator haploCount= input1.begin();haploCount !=input1.end();++haploCount)
   {
      vector<std::pair<int,int>> my_map; // new sub container      
      std::transform((*haploCount).begin(), (*haploCount).end(), input2.begin(), 
      std::inserter(my_map, my_map.end()), std::make_pair<int const&,int const&>);// push in a map for a single haplo vector and the state vector
		      
      std::sort(my_map.begin(),my_map.end(),pairCompare); // sort the map      
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