#include "particlefilter.h"


particle_filter::particle_filter()
{
  
}
particle_filter::~particle_filter()
{
  
}


void particle_filter::filter_run(const tbb::concurrent_vector<vector<vector<Double>>> &Forward_values, 
		       const tbb::concurrent_vector<vector<vector<Double>>> &Backward_values,
		       const vector<vector<vector<int>>> &curStateSpace,const vector< states > &tot_StateSpace,
		       vector<vector<vector<int>>> &upd_StateSpace)
{
  //cout << "particle_filter::filter_run Start" << endl;
  
  vector<vector<Double>> sum_ForwdVals, sum_BckwdVals;
  
  vector<vector<int>> states_sortedFwd,states_sortedBwd,states_Current;
  
  //first fetch the summary of all individuals for forward and the reverse probabilities
  //summary_Markervals(Forward_values,sum_ForwdVals);
  //summary_Markervals(Backward_values,sum_BckwdVals);
  
  for(int indCount=0; indCount<n_individuals;++indCount)
  {
      states_Current = curStateSpace[indCount];
      vector<vector<int>> indData;
      // Sort the states at all markers based on Fwd and back probs
      sum_ForwdVals = Forward_values[indCount];
      sum_BckwdVals = Backward_values[indCount];
      rank_Statescontrib(sum_ForwdVals,states_Current,states_sortedFwd);
      rank_Statescontrib(sum_BckwdVals,states_Current,states_sortedBwd);
  
      //Finally filter out the poor states and include prospective states
      filter_states(states_sortedFwd,states_sortedBwd,tot_StateSpace,indData);
      upd_StateSpace.emplace_back(indData);
      states_sortedFwd.clear();
      states_sortedBwd.clear();
  }
  //cout << updStateSpace.size() << endl;
 // cout << "particle_filter::filter_run End" << endl;
}


void particle_filter::summary_Markervals(const vector<vector<vector<Double>>> &Input,vector<vector<Double>> &T)
{
  Double stateTotal;
   
   for(int markCount = 0;markCount < n_markers;markCount++)   
   {
	vector<Double > markerData;
	
	for(int stateCount=0;stateCount<num_states;stateCount++)
	{
	    stateTotal = 0.0;
	    for(int indCount = 0;indCount < n_individuals;indCount++)
	    {
		stateTotal = stateTotal + Input[indCount][markCount][stateCount];
	    }
	    markerData.push_back(stateTotal);
     	}
     	
	T.push_back(markerData);
   }
   
}

void particle_filter::rank_Statescontrib(vector<vector<Double>> &Input,
	const vector<vector<int>> &curStateSpace, vector<vector<int>> &T)
{
  
 // cout << "particle_filter::rank_Statescontrib Start" << endl;
  
  //ofstream myfile;
  //myfile.open ("Phase_p2.txt");  
  
  vector< multimap<Double,int>> mappedVals;  
  vector<vector<Double>>::iterator ValsCount= Input.begin();
  
  //first pair the state and the value at the state in a map container to obtain a sorted map in ascending order
  for(vector<vector<int>>::const_iterator markCount= curStateSpace.begin();markCount !=curStateSpace.end();markCount++,ValsCount++)
  {
       multimap<Double,int> my_map;
       std::transform((*ValsCount).begin(), (*ValsCount).end(), (*markCount).cbegin(), 
       std::inserter(my_map, my_map.end()), std::make_pair<Double const&,int const&>);
       
      mappedVals.emplace_back(my_map);
  }
 
   
  /*for(multimap<Double, int>::iterator it = mappedVals[0].begin(); it != mappedVals[0].end(); it++)
        cout << it -> first << " " << it -> second << endl; */
  
  //copy the sorted states into the output variables
  for(vector<multimap<Double,int>>::iterator markCount = mappedVals.begin(); markCount != mappedVals.end();markCount++)
  {
      vector<int> markerData;
        
      for(multimap<Double,int>::const_reverse_iterator it = markCount->rbegin(); it != markCount->rend(); ++it)
      {
	  //myfile << it -> first << " " << it->second << endl; 
	  markerData.emplace_back(it -> second);
      }
      
      //myfile << endl<< " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
      T.emplace_back(markerData);
  }
    
   //myfile.close(); 
 //  cout << "particle_filter::rank_Statescontrib End" << endl;
}


void particle_filter::filter_states(vector<vector<int>> &sortedStsFwd,vector<vector<int>> &sortedStsBwd,
			 const vector<states>tot_Space,vector<vector<int>> &upd_StateSpace)
{
  
    //cout << "particle_filter::filter_states Start" << endl;
    int iStateskeep ;    
   
    vector<int> statesKeep;      
    
    iStateskeep = (sortedStsFwd[0].size()/3); //keep one third the number of states
   
    int i=0,iCounter,totalstates;
    vector<states > totSpaceCopy;
    totalstates= tot_Space.size();
    vector<int>::iterator marker_it;   
    vector<int>::const_reverse_iterator rev_marker_it ;
    
    for(int iMarkCount=0;iMarkCount< n_markers;++iMarkCount)
    {
      totSpaceCopy = tot_Space;
      
      //take the best  states with fwdvalues
      statesKeep = vector<int>{sortedStsFwd[iMarkCount].begin(),sortedStsFwd[iMarkCount].begin()+ iStateskeep};
     
      //take the best  states with bckwrdvalues
       marker_it = sortedStsBwd[iMarkCount].begin();
       
      for(int iCounter=0;iCounter < iStateskeep ;++iCounter)
      {
	if(std::find(statesKeep. begin(),statesKeep.end(),(*marker_it))== statesKeep.end())
	  statesKeep.emplace_back(*marker_it);
	++marker_it;
      }
      vector<int> markData{statesKeep};
      
      statesKeep.clear();
      //elimininate the states with bad fwdvalues
      statesKeep.insert(statesKeep.end(),sortedStsFwd[iMarkCount].rbegin(),sortedStsFwd[iMarkCount].rbegin()+iStateskeep);
      
      //elimininate the states with bad bckwrdvalues
      rev_marker_it = sortedStsBwd[iMarkCount].rbegin();
      for(int iCounter=0;iCounter < iStateskeep ;++iCounter)
      {
	if(std::find(statesKeep.begin(),statesKeep.end(),(*rev_marker_it))== statesKeep.end())
	  statesKeep.emplace_back(*rev_marker_it);
	
	++rev_marker_it;
      }
      std::sort(statesKeep.begin(),statesKeep.end());
      statesKeep.erase( unique( statesKeep.begin(), statesKeep.end() ), statesKeep.end());              
      sample(totSpaceCopy,statesKeep,markData,num_states);
      statesKeep.clear();
      totSpaceCopy.clear();
      upd_StateSpace.push_back(markData);       
    }
     
    /*ofstream myfile;
    myfile.open ("Phase_filter.txt",ios::app|ios::in);   
    myfile << "sortedStsFwd" << endl;
    
    for(vector<vector<int>>::iterator markCountf= sortedStsFwd.begin();markCountf !=sortedStsFwd.end();markCountf++)
    {
	for(vector<int>::iterator stateCountf = markCountf->begin();stateCountf != markCountf->end();stateCountf++)
	{
	  myfile << (*stateCountf) << " "; 
	}
	myfile << endl;
    }
     
   myfile << "sortedStsBwd" << endl;
   for(vector<vector<int>>::iterator markCount= sortedStsBwd.begin();markCount !=sortedStsBwd.end();markCount++)
   {
      for(vector<int>::iterator stateCount= markCount->begin();stateCount != markCount->end();stateCount++)
      {
	  myfile << (*stateCount) << " "; 
      }
      myfile << endl;
    }
    
    myfile << "updated states" << endl;
    for(vector<vector<int>>::iterator markCountf= upd_StateSpace.begin();markCountf !=upd_StateSpace.end();markCountf++)
    {
	for(vector<int>::iterator stateCountf = markCountf->begin();stateCountf != markCountf->end();stateCountf++)
	{
	  myfile << (*stateCountf) << " "; 
	}
	myfile << endl;
    }
   myfile.close(); */
  // cout << "particle_filter::filter_states End" << endl;
   
}


















