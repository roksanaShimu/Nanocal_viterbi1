#ifndef __PARALLEL_IMPLEMENTATION_1i_HPP
#define __PARALLEL_IMPLEMENTATION_1i_HPP


//m  read: texture memory
//m  write: global memory
// n_ from v: page lock memory
// pm : page lock memory
//ev: constant memory 


#include<iostream>
#include<vector>
#include<stdlib.h>
#include <utility>
#include <math.h>
#include<stdio.h>

#include "Data_Structure.hpp"
#include "Predefine_Values.hpp"
#include "Device_Memory_Allocation.hpp"

using namespace std;


texture<float, cudaTextureType1D, cudaReadModeElementType> texMFloat;
texture<unsigned, cudaTextureType1D, cudaReadModeElementType> texMUnsigned;
__constant__ float d_ev[4];


__global__ void kernel_calculation(unsigned * m_curr_unsigned, float * m_curr_float, unsigned * devPtr_n_fv_unsigned, float * devPtr_n_fv_float, float * devPtr_pm, unsigned n_states){

	//texMFloat: m previous event float
	//texMUnsigned: m previous event float


	unsigned id = threadIdx.x + blockIdx.x * threads_per_block;   //id is j from the serial code

	
	float m_alpha=-INFINITY;
	unsigned m_beta=n_states;
	
	unsigned group_no=id/warp_size;
	__syncthreads();
	  
	for(unsigned k=0;k<length_of_from_v;k++){

		const unsigned& j_prev = devPtr_n_fv_unsigned[group_no*warp_size*length_of_from_v + warp_size*k + id%warp_size];
                const float& log_pr_transition = devPtr_n_fv_float[group_no*warp_size*length_of_from_v + warp_size*k + id%warp_size];
		
                float v = log_pr_transition + tex1D(texMFloat, j_prev);

                if (v > m_alpha)
                    {
                        m_alpha = v;
                        m_beta = j_prev;
                    }
                
        }

	__syncthreads();
	unsigned pm_index= group_no* (warp_size*4)+ (id%warp_size);  //index of first element

	float log_level_stdv=log( devPtr_pm[pm_index+warp_size]);     

	float log_2pi = log(2.0 * M_PI);
	float a=(d_ev[0] - devPtr_pm[pm_index] ) / devPtr_pm[pm_index+warp_size];
	float log_normal_pdf=- log_level_stdv - (log_2pi + a * a) /static_cast< float >(2.0);


	float log_stdv=log(d_ev[1]);
	float sd_lambda = pow(devPtr_pm[pm_index+2*warp_size], 3.0f) / pow(devPtr_pm[pm_index+3*warp_size], 2.0f);
	float log_sd_lambda=log(sd_lambda );     // do i need to scale them?

	
	float b = (d_ev[1] - devPtr_pm[pm_index+2*warp_size]) / devPtr_pm[pm_index+2*warp_size];
	float log_invgauss_pdf=(log_sd_lambda - log_2pi - static_cast< float >(3.0) * log_stdv - sd_lambda * b * b / d_ev[1]) / static_cast< float >(2.0);

	m_alpha += log_normal_pdf+log_invgauss_pdf;

	__syncthreads();
	m_curr_unsigned[id]=m_beta;
	m_curr_float[id]=m_alpha;

}


void update_the_m(unsigned * h_m_curr_unsigned, float * h_m_curr_float, vector<Matrix_Entry> & m, unsigned event_number, unsigned n_states){
	unsigned start_index=event_number*n_states; 
	for(unsigned j=0; j<n_states; j++){
		m[start_index+j].alpha=h_m_curr_float[j];
		m[start_index+j].beta=h_m_curr_unsigned[j];
	}
}


void parallel_calculations_1(unsigned n_states, vector<Matrix_Entry> & m , vector< vector<pair<unsigned, float> > > neighbors_from_v, vector<event> ev, vector<pore_model_state> pm){

	// n_states is 4096
	// n_events= ev.size() // = 2972
	// m : size= n_events * n_states   //for 1 event size= n_states * 8 bytes = 32768 B // it also requires info of previous event
	// neighbors_from_v = n_states * length_of_from_v * 8 Bytes   //4096 * 21 * 8B= 32768 B * 21 = 688128  //exceeds all memories
	// ev: n_events * 16 Bytes
	// pm: n_states* 16 Bytes  // 4096* 16B= 65536 B //reaches the limit of most of the memory limit

	
	//experiment 1:
	/*	for i= each of the events {  //2972 times
			
			keep m(i-1)(n_states) in any read_only_memory and event dependent  // or make it surface memory 
			allocate space in the global memory for m(i)(n_states). results will be written here by the threads
			neighbors_from_v is the biggest read only data. So use page lock option for this variable //remains same for all event
			pm also read_only data and event independant. keep it in constant memory //though it will require the whole memory
			ev

		}
	*/



	//*****************************initialize device parameter**********************************

	int devId = 0;
	//if (argc > 1) devId = atoi(argv[1]);

	cudaDeviceProp prop;
	checkCuda( cudaGetDeviceProperties(&prop, devId));
	printf("Device : %s\n", prop.name);
	checkCuda( cudaSetDevice(devId) );

	cudaGetDeviceProperties(&prop, 0);
	if (!prop.canMapHostMemory) {
		printf("Device %d cannot map host memory!\n", 0);
		exit(EXIT_FAILURE);

	}/*else{
		cout<<"can be mapped"<<endl;
	}*/

	//***************************************************************

	//allocate texture read_only memory for m(i-1)(n_states): m_pre
	//n_states thread will read n_states data. So there is no need to rearrange the data

	float* m_float; unsigned * m_unsigned;
	
	/* Allocate space for results on host */
    	m_float = (float *)valloc(n_states*sizeof(float));
	m_unsigned = (unsigned *)valloc(n_states*sizeof(unsigned));
	

	for(unsigned i=0;i<n_states;i++){ // for event 0
		m_float[i]=m[i].alpha;
		m_unsigned[i]=m[i].beta;
	}	


	allocate_texture_memory_float(texMFloat, m_float, n_states*sizeof(float), n_states); 
	allocate_texture_memory_unsigned(texMUnsigned, m_unsigned, n_states*sizeof(unsigned), n_states); 
//--------------------------------------------
	

//********************************************************************************


	// allocate global memory for m(i)(n_states): m_curr
	unsigned  *m_curr_unsigned; float  *m_curr_float;
	cudaMalloc((void**)&m_curr_unsigned, sizeof(unsigned)*n_states) ; // device
	cudaMalloc((void**)&m_curr_float, sizeof(float)*n_states) ; // device

	//unsigned  *h_m_curr_unsigned; float  *h_m_curr_float;
	//h_m_curr_unsigned = (unsigned *)valloc(sizeof(unsigned) * n_states);
	//h_m_curr_float = (float *)valloc(sizeof(float) * n_states);
	

	
//********************************************************************************

	//page lock memory for neighbors_from_v //event independant  :n_from_v
	//rearrange the values in such a way that memory reading will be coalesced
	// neighbors_from_v have n_states and each state has "length_of_from_v" elements
	// Initially the 2d vector is organized as follow:
	//	state[0].from_v[0] 	state[0].from_v[1]	state[0].from_v[2]  ...  state[0].from_v[length_of_from_v-1]
	//	state[1].from_v[0]	state[1].from_v[1]	state[1].from_v[2]  ...	 state[1].from_v[length_of_from_v-1]
	//		.			.			.	    ...		.
	//		.			.			.           ...		.
	// The vector is reorganized to a 1D array as:
	//	state[0].from_v[0] 	state[1].from_v[0]	state[2].from_v[0]  ...  state[warp.size-1].from_v[0] 	state[0].from_v[1]	state[1].from_v[1]	state[2].from_v[1]  ...	 state[warp.size-1].from_v[1] .....

	unsigned * n_fv_unsigned; float * n_fv_float;
	size_t n_fv_unsigned_size= sizeof(unsigned) * n_states * length_of_from_v;
	size_t n_fv_float_size= sizeof(float) * n_states * length_of_from_v;

	/* Allocate space on host */
    	n_fv_unsigned = (unsigned *)valloc(n_fv_unsigned_size);
	n_fv_float = (float *)valloc(n_fv_float_size);
	
	for(unsigned i=0; i<neighbors_from_v.size(); i++){
		for(unsigned j=0;j<neighbors_from_v[i].size();j++){
			unsigned group_no=i/warp_size;
			unsigned index= group_no* (warp_size*length_of_from_v)+ warp_size*j + (i%warp_size);
			n_fv_unsigned[index]=neighbors_from_v[i][j].first; 
			n_fv_float[index]=neighbors_from_v[i][j].second; 

		}
	}

	; 
	unsigned *devPtr_n_fv_unsigned=assign_page_locked_memory_unsigned( n_fv_unsigned, n_fv_unsigned_size);

	float *devPtr_n_fv_float=assign_page_locked_memory_float(n_fv_float, n_fv_float_size);
	
	

//********************************************************************************

	//allocate page_lock memory for pm //event independant   : pm

	// arrange the values so that the memory access time will be the least
	float * pm_rearranged; size_t pm_rearranged_size=n_states * 4* sizeof(float);
	pm_rearranged = (float *)valloc(pm_rearranged_size);

	for(unsigned i=0; i<n_states; i++){
		unsigned group_no=i/warp_size;
		unsigned index= group_no* (warp_size*4)+ (i%warp_size);  //index of first element
		unsigned n=warp_size; //theads_belong_to_this_group                          //do the logics
 		pm_rearranged[index]=pm[i].level_mean; pm_rearranged[index+n]= pm[i].level_stdv;
		pm_rearranged[index+2*n]=pm[i].sd_mean; pm_rearranged[index+3*n]= pm[i].sd_stdv;
	} 
	

	float *devPtr_pm= assign_page_locked_memory_float(pm_rearranged, pm_rearranged_size);
	

	cout<<"allocations are completed"<<endl;

//********************************************************************************

	unsigned n_events=ev.size();
	//cout<<"n_events = "<<n_events<<endl;
	for(unsigned i=1; i<n_events; i++){

		//assign ev in constant memory
		float ev_arr[4]; ev_arr[0]=ev[i].mean;  ev_arr[1]=ev[i].stdv;  ev_arr[2]=ev[i].start;  ev_arr[3]=ev[i].length;  
		cudaMemcpyToSymbol(d_ev, ev_arr, 4*sizeof(float));

		//assign m_pre
		allocate_texture_memory_float(texMFloat, m_float, n_states*sizeof(float), n_states); 
		//allocate_texture_memory_unsigned(texMUnsigned, m_unsigned, n_states*sizeof(unsigned), n_states);

		//clock_t kernel_start = clock();
		//call kernel (m_pre, m_curr, n_from_v, pm, ev)
		kernel_calculation<<< number_of_blocks, threads_per_block>>>(m_curr_unsigned, m_curr_float, devPtr_n_fv_unsigned, devPtr_n_fv_float, devPtr_pm, n_states );
		cudaDeviceSynchronize();
		
		//copy m_cuur	
		cudaMemcpy(m_unsigned,m_curr_unsigned, sizeof(unsigned)*n_states, cudaMemcpyDeviceToHost);
		cudaMemcpy(m_float,m_curr_float, sizeof(float)*n_states, cudaMemcpyDeviceToHost);
		
		cudaUnbindTexture(texMFloat);

		update_the_m(m_unsigned, m_float, m, i, n_states);
		//clock_t kernel_end = clock();
		//printf("Time taken for kernel calculation: %.6fs\n", (double)(kernel_end - kernel_start)/CLOCKS_PER_SEC);

		//cout<<"i="<<i<<endl;
	}

	


}




#endif
