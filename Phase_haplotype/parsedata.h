#pragma once
#define PARSEDATA_H_INCLUDED

#include "globals.h"




class data_genotypes
{
    public:
    
    vector<vector<vector<int>>>   m_genotypes; //nIndividuals x nmarkers x nhaplotypes
    vector<double> m_physical_positions;
    vector<double> m_physical_distances;
    //vector<std::string> m_ind_ID;
    void read_paramsFile();
    data_genotypes();
    ~data_genotypes();
    void parse_input(string ipFileName);
    void print_variables(void);
};