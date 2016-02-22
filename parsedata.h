#pragma once
#define PARSEDATA_H_INCLUDED

#include "globals.h"




class data_genotypes
{
    public:
    
    vector<vector<vector<int>>>   m_genotypes; //nIndividuals x nmarkers x nhaplotypes
    vector<int> m_physical_positions;
    vector<int> m_physical_distances;
    vector<std::string> m_ind_ID;
    
    data_genotypes();
    ~data_genotypes();
    void parse_input(void);
    void log_variables(void);
    void read_variables();
};