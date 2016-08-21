

#include<iostream>
#include<vector>
#include<stdlib.h>
#include <utility>

#include "Initialization.hpp"
#include "Predefine_Values.hpp"
#include "Data_Structure.hpp"
#include "Serial_Implementation.hpp"
#include "Parallel_Implementation_1i.hpp"
//#include "Parallel_Implementation_1ii.hpp"
//#include "Parallel_Implementation_1ib.hpp"

using namespace std;


int main(){
	unsigned n_states=4096*2;
	unsigned n_events=2972*1;
	vector<Matrix_Entry> m; m.resize(n_states*n_events);  // need to calculate the values of m
	//fill m for first event;
	fill_randomly_m_for_first_event(m, n_states, n_events);

	//for error checking 
	vector<Matrix_Entry> m1; m1.resize(n_states*n_events);m1=m;
	vector<Matrix_Entry> m2; m2.resize(n_states*n_events);m2=m;

	vector< vector<pair<unsigned, float> > > neighbors_from_v; neighbors_from_v.resize(n_states);
	for(unsigned i=0; i<n_states; i++){
		neighbors_from_v[i].resize(length_of_from_v);
	}
	generate_random_values_for_neighbors_from_v(neighbors_from_v);
	//print_neighbors_from_v(neighbors_from_v);
	

	vector<pore_model_state> pm; pm.resize(n_states);
	generate_random_values_for_pm(pm, n_states);

	
	vector<event> ev; ev.resize(n_events);
	generate_random_values_for_event(ev);
	//print_event(ev);
	

	clock_t serial_start = clock();
	serial_calculations(n_states, m1, neighbors_from_v, ev, pm);
	clock_t serial_end = clock();

	cout<< "serial calculation is over"<<endl;
	printf("Time taken for serial_code: %.6fs\n", (double)(serial_end - serial_start)/CLOCKS_PER_SEC);

	clock_t parallel_start = clock();
	parallel_calculations_1(n_states, m2, neighbors_from_v, ev, pm);
	clock_t parallel_end = clock();


	clock_t parallel_1ii_start = clock();
	//parallel_calculations_1ii(n_states, m2, neighbors_from_v, ev, pm);
	clock_t parallel_1ii_end = clock();

	clock_t parallel_1ib_start = clock();
	//parallel_calculations_1ib(n_states, m2, neighbors_from_v, ev, pm);
	clock_t parallel_1ib_end = clock();


	printf("Time taken for serial_code: %.6fs\n", (double)(serial_end - serial_start)/CLOCKS_PER_SEC);
	printf("Time taken for parallel: %.6fs\n", (double)(parallel_end - parallel_start)/CLOCKS_PER_SEC);
	printf("Time taken for parallel: %.6fs\n", (double)(parallel_1ii_end - parallel_1ii_start)/CLOCKS_PER_SEC);
	printf("Time taken for parallel: %.6fs\n", (double)(parallel_1ib_end - parallel_1ib_start)/CLOCKS_PER_SEC);



 
	int error=0;
	int counter=0;
	for(unsigned i=0; i<m.size(); i++){
		if(m1[i].beta != m2[i].beta){
			//cout<<i<<"th value didn't matched. m1[i].beta= "<<m1[i].beta<<",  m1[i].alpha= "<< m1[i].alpha<< "  and m2[i].beta= "<< m2[i].beta<<", m2[i].alpha= "<<m2[i].alpha<<endl;
			error=1;
			counter++;
			//break;
		}
	}
	if(error==0){
		cout<<"all matched"<<endl;
	}else{
		cout<<counter<<" elements didn't matched"<<endl;
	}
	return 0;
}

