#include "particlefilter.h"


particle_filter::particle_filter()
{
  
}
particle_filter::~particle_filter()
{
  
}

//############################## WHEN individuals share same state space and states are same for all markers###############
//#################################################################################
void  particle_filter::filter_run(const vector<vector<vector<Double>>> &hmm_values,const vector<int> &curStateSpace,
				  const vector< Double > &tot_StateSpace_weights,  vector<int> &upd_StateSpace,Chaplotypes &HmmObj)
{
    //cout << "particle_filter::filter_run Start" << endl;
    
    vector<Double> sum_AggVals;  
    vector<int> states_sorted;
    
    //first fetch the summary of all individuals for forward and the reverse probabilities
    summary_Markervals(hmm_values, sum_AggVals);
  
    // Sort the states at all markers based on Fwd and back probs
    rank_Statescontrib(sum_AggVals,curStateSpace,states_sorted);
    filter_states(states_sorted,tot_StateSpace_weights,upd_StateSpace,HmmObj);

   //cout << "particle_filter::filter_run End" << endl;
}
 void particle_filter::summary_Markervals(const vector<vector<vector<Double>>> &hmm_values,  vector<Double> &T)
{
    
      vector<vector<Double>>  allMarkerData;
      T= vector<Double>(num_states);
      for(int markCount = 0;markCount < n_markers;++markCount)   
      {
	    vector<Double> markerData(num_states);
	  
	     parallel_for(int(0), num_states, [this,&markerData,&hmm_values,&markCount](int stateCount)throw()
	    {
		  Double stateTotal = 1.0;
		  for(int indCount = 0;indCount < n_individuals;++indCount)
		  {
		      stateTotal *= hmm_values[indCount][markCount][stateCount];
		  }
		  markerData[stateCount] = stateTotal;		
	    });
	    
	   allMarkerData.push_back(markerData);
      }
       
	parallel_for(int(0), num_states, [this,&allMarkerData,&T](int stateCount)throw()
	{
	      Double stateTotal = 0.0;
	      for(int markCount = 0;markCount < n_markers;++markCount) 
	      {
		    stateTotal += allMarkerData[markCount][stateCount];	    
	      }
	      T[stateCount]= stateTotal;
      });
}

//retains the top one third states  and samples for new states surrounding these states after discarding states with low values
void particle_filter::filter_states(vector<int> &sortedSts, const vector<Double>&tot_Space,vector<int> &upd_StateSpace,Chaplotypes &HmmObj)
{
      //cout << "particle_filter::filter_states Start" << endl;
      int iStateskeep = (num_states/3); //keep one third the number of states      
      int inew_states = num_states - (iStateskeep);
      int itotal_states= tot_Space.size();
      std::vector<int> upd_StateSpace_temp;
      
      upd_StateSpace =std::vector<int>(iStateskeep);      
      //ADD good states which we want to retain in next iteration
      std::transform( sortedSts.begin(), sortedSts.begin()+(iStateskeep), upd_StateSpace.begin(),
	      [] (int const& ms){ return ms;}); 
    
      std::vector<int> statesExclude(iStateskeep*2);      
      //save good states which we want to retain in next iteration, to avoid re-sampling them at the end
      std::transform( sortedSts.begin(), sortedSts.begin()+(iStateskeep), statesExclude.begin(),
	      [] (int const& ms){ return ms;});     
      //save bad states which we want to eliminate from next iteration, to avoid re-sampling them at the end
      std::transform( sortedSts.rbegin(), sortedSts.rbegin()+(iStateskeep), statesExclude.begin()+(iStateskeep),
	      [] (int const& ms){ return ms;});
     
      std::sort(statesExclude.begin(),statesExclude.end());
   
      //keep track of remaining states after excluding good and bad states
      std::vector<int> remaining_states (itotal_states);
      std::generate_n (remaining_states.begin (),itotal_states, [] { static int i {0}; return i++; });    
    
     remaining_states.erase( std::remove_if(remaining_states.begin(), remaining_states.end(), 
	[&statesExclude](int i){ return (std::find(statesExclude.begin(), statesExclude.end(), i) != statesExclude.end());}),remaining_states.end());
           
      std::sort(statesExclude.begin(),statesExclude.end(), std::greater<int>());//sort in descending order      
      // exclude the elements from total state space before sampling for new states for next iteration
      std::vector<Double > totSpaceCopy = tot_Space; 
      for ( int i : statesExclude ) totSpaceCopy.erase( std::next( totSpaceCopy.begin(), i ) );
   
      sample(totSpaceCopy,upd_StateSpace_temp,inew_states);    
      
      //retrieve the original indices with out the impact of excluded states
      for(auto &it:upd_StateSpace_temp)
      {
	  upd_StateSpace.emplace_back(remaining_states[it]);
      }
     // cout << "particle_filter::filter_states End" << endl;
}

 //ranks the states sorted according the summary of values for all individuals at each markerstate
void particle_filter::rank_Statescontrib(vector<Double> &Input, const vector<int> &curStateSpace, vector<int> &T)
{
    // cout << "particle_filter::rank_Statescontrib Start" << endl;
  
   multimap<Double,int> mappedVals;  
  
  //first pair the state and the value at the state in a map container to obtain a sorted map in ascending order
    std::transform(Input.begin(), Input.end(), curStateSpace.cbegin(), 
    std::inserter(mappedVals, mappedVals.end()), std::make_pair<Double const&,int const&>);
 
  //copy the sorted states into the output variables
  for(multimap<Double,int>::const_reverse_iterator it = mappedVals.rbegin(); it != mappedVals.rend(); ++it)
  {
	T.emplace_back(it -> second);
  }
 
 //  cout << "particle_filter::rank_Statescontrib End" << endl;
}




