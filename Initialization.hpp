#ifndef __INITIALIZATION_HPP
#define __INITIALIZATION_HPP

#include<iostream>
#include<vector>
#include<stdlib.h>
#include <utility>

#include "Data_Structure.hpp"

using namespace std;


//*************************************************************
void fill_randomly_m_for_first_event(vector<Matrix_Entry>&m, unsigned n_states, unsigned n_events){
	if(n_events==0)return;
	for(unsigned i=0;i<n_states;i++){
		m[i].alpha=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		m[i].beta= (unsigned) rand()%4096;
	}
}

//*************************************************************
void generate_random_values_for_neighbors_from_v(vector< vector<pair<unsigned, float> > > & neighbors_from_v){
	typedef pair<unsigned, float> intpair;
	for(int i=0; i<neighbors_from_v.size(); i++){
		for(unsigned j=0;j<neighbors_from_v[i].size();j++){
			float r= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			unsigned n= (unsigned) rand()%4096;
			neighbors_from_v[i][j]=intpair(n,r);
		}
		
	}
}
void print_neighbors_from_v(vector< vector<pair<unsigned, float> > > neighbors_from_v){
	for(unsigned i=0;i<neighbors_from_v.size();i++){
		for(unsigned j=0; j<neighbors_from_v[i].size(); j++){
			cout<<neighbors_from_v[i][j].first<<"   "<<neighbors_from_v[i][j].second<<endl;
		}
		cout<<endl;
	}
}

//*************************************************************
void generate_random_values_for_pm(vector<pore_model_state>& pm, unsigned n_states){
	for(unsigned i=0; i<n_states; i++){
		
		pm[i].level_mean=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		pm[i].level_stdv=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		pm[i].sd_mean=static_cast <float> (rand()) / static_cast <float> (RAND_MAX); 
		pm[i].sd_stdv=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	}
}
//*************************************************************
void generate_random_values_for_event(vector<event> & ev){
	for(int i=0; i<ev.size(); i++){
		ev[i].mean=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		ev[i].stdv=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		ev[i].start= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		ev[i].length=static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		
	}
}
void print_event(vector<event> ev){
	for(int i=0; i<ev.size(); i++){
		cout<<ev[i].mean<<"  "<<ev[i].stdv<<"   "<<ev[i].start<<"   "<<ev[i].length<<endl;
	}

}




#endif
