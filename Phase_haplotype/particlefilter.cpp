#include "particlefilter.h"


particle_filter::particle_filter()
{
  
}
particle_filter::~particle_filter()
{
  
}


void particle_filter::filter_run(const vector<vector<vector<Double>>> &hmm_values, 
		        const vector<vector<vector<int>>> &curStateSpace,const vector< Double > &tot_StateSpace_weights,
		       vector<vector<vector<int>>> &upd_StateSpace)
{
  //cout << "particle_filter::filter_run Start" << endl;
  
	vector<vector<Double>> ind_hmmVals;
	
	vector<vector<int>> states_sortedhmm,states_Current;
	
	//first fetch the summary of all individuals for forward and the reverse probabilities
		
	for(int indCount=0; indCount<n_individuals;++indCount)
	{
		  states_Current = curStateSpace[indCount];
		  ind_hmmVals = hmm_values[indCount];
		  vector<vector<int>> indData;
		  // Sort the states at all markers based on Fwd and back probs
		  rank_Statescontrib(ind_hmmVals,states_Current,states_sortedhmm);
		  	      
		  //Finally filter out the poor states and include prospective states
		  filter_states(states_sortedhmm,tot_StateSpace_weights,indData);
		  upd_StateSpace.emplace_back(indData);
		  states_sortedhmm.clear();
	}
 // cout << "particle_filter::filter_run End" << endl;
}



void particle_filter::rank_Statescontrib(vector<vector<Double>> &Input,
	const vector<vector<int>> &curStateSpace, vector<vector<int>> &T)
{
  
 // cout << "particle_filter::rank_Statescontrib Start" << endl;
      multimap<Double,int> mappedVals;  
      for(int markCount=0;markCount<n_markers;++markCount)
      {
	      vector<int> marker_values;
	      std::transform(Input[markCount].begin(), Input[markCount].end(), curStateSpace[markCount].cbegin(), 
	      std::inserter(mappedVals, mappedVals.end()), std::make_pair<Double const&,int const&>);
	    
	      //copy the sorted states into the output variables
	      for(multimap<Double,int>::const_reverse_iterator it = mappedVals.rbegin(); it != mappedVals.rend(); ++it)
	      {
		   marker_values.emplace_back(it -> second);
	      }
	      T.emplace_back(marker_values);
	      mappedVals.clear();
      }
      

 //  cout << "particle_filter::rank_Statescontrib End" << endl;
}

void particle_filter::filter_states(vector<vector<int>> &sortedSts,
			 const vector<Double>tot_Space_weights,vector<vector<int>> &upd_StateSpace)
{
  
    //cout << "particle_filter::filter_states Start" << endl;
	int iStateskeep = (num_states/3); //keep one third the number of states      
	int inew_states = num_states - (iStateskeep);
	int itotal_states= tot_Space_weights.size();
	std::vector<int> upd_StateSpace_temp;
	std::vector<int> statesExclude;
	std::vector<int> remaining_states ;
	std::vector<Double > totSpaceCopy;
	vector <int> state_ids(itotal_states); 
	int stateInd=-1;
	for(auto &it:state_ids)
	{
	      it = ++stateInd;
	}
      for(int markCount=0;markCount< n_markers;++markCount)
     {
	     vector<int> upd_markerStates(iStateskeep); 
	   
	    //ADD good states which we want to retain in next iteration
	      std::transform( sortedSts[markCount].begin(), sortedSts[markCount].begin()+(iStateskeep), upd_markerStates.begin(),
	      [] (int const& ms){ return ms;}); 
	    
	      statesExclude = std::vector<int>(iStateskeep*2);
	    
	     //save good states which we want to retain in next iteration, to avoid re-sampling them at the end
	      std::transform(  sortedSts[markCount].begin(),  sortedSts[markCount].begin()+(iStateskeep), statesExclude.begin(),
	      [] (int const& ms){ return ms;});   
	      
	      //save bad states which we want to eliminate from next iteration, to avoid re-sampling them at the end
	      std::transform( sortedSts[markCount].rbegin(), sortedSts[markCount].rbegin()+(iStateskeep), statesExclude.begin()+(iStateskeep),
	      [] (int const& ms){ return ms;});
	      std::sort(statesExclude.begin(),statesExclude.end());
	      
	        //keep track of remaining states after excluding good and bad states
	      remaining_states =state_ids;
	
	       remaining_states.erase( std::remove_if(remaining_states.begin(), remaining_states.end(), 
	      [&statesExclude](int i){ return (std::find(statesExclude.begin(), statesExclude.end(), i) != statesExclude.end());}),remaining_states.end());
	   
	       std::sort(statesExclude.begin(),statesExclude.end(), std::greater<int>());//sort in descending order
	       
	         // exclude the elements from total state space before sampling for new states for next iteration
		totSpaceCopy = tot_Space_weights; 
		for ( int i : statesExclude ) totSpaceCopy.erase( std::next( totSpaceCopy.begin(), i ) );
		
		sample(totSpaceCopy,upd_StateSpace_temp,inew_states);   
		
		 //retrieve the original indices with out the impact of excluded states
		for(auto &it:upd_StateSpace_temp)
		{
		     upd_markerStates.emplace_back(remaining_states[it]);
		}
	
		upd_StateSpace.emplace_back(upd_markerStates);
		upd_StateSpace_temp.clear();
		remaining_states.clear();
    }
   
}
















