#ifndef __SERIAL_IMPLEMENTATION_HPP
#define __SERIAL_IMPLEMENTATION_HPP


#include<iostream>
#include<vector>
#include<stdlib.h>
#include <utility>

#include "Data_Structure.hpp"

using namespace std;

void serial_calculations(unsigned n_states, vector<Matrix_Entry> & m , vector< vector<pair<unsigned, float> > > neighbors_from_v, vector<event> ev, vector<pore_model_state> pm){
	unsigned n_events=ev.size();
	for (unsigned i = 1; i < n_events; ++i)
        {
            
            for (unsigned j = 0; j < n_states; ++j) // TODO: parallelize
            {
                m[i*n_states+j].alpha = -INFINITY;
                m[i*n_states+j].beta = n_states;
               	
                for(unsigned k=0;k<neighbors_from_v[j].size();k++){

                    const unsigned& j_prev = neighbors_from_v[j][k].first;
                    const float& log_pr_transition = neighbors_from_v[j][k].second;

                    float v = log_pr_transition + m[(i - 1)*n_states+ j_prev].alpha;

                    if (v > m[(i*n_states)+ j].alpha)
                    {
                        m[(i*n_states)+ j].alpha = v;
                        m[(i*n_states)+ j].beta = j_prev;
                    }
                
                }
		/*m[i*n_states+j].alpha += pm.log_pr_emission(j, ev[i]);*/
	
                // for log_normal pdf	
		/*float ev[i].mean;	
		float level_mean=???
		float level_stdv=?*/
		float log_level_stdv=log(pm[j].level_stdv);

		static const float log_2pi = log(2.0 * M_PI);
		float a=(ev[i].mean - pm[j].level_mean) / pm[j].level_stdv;
		float log_normal_pdf=- log_level_stdv - (log_2pi + a * a) / static_cast< float >(2.0);


		//for log_invgauss_pdf
		//float ev[i].stdv;
		//float sd_mean=ev[i].sd_mean??
		//float sd_stdv=ev[i].sd_stdv?

		float log_stdv=log(ev[i].stdv);
		float sd_lambda = pow(pm[j].sd_mean, 3.0) / pow(pm[j].sd_stdv, 2.0);
		float log_sd_lambda=log(sd_lambda );     // do i need to scale them?
		
		float b = (ev[i].stdv - pm[j].sd_mean) / pm[j].sd_mean;
		float log_invgauss_pdf=(log_sd_lambda - log_2pi - static_cast< float >(3.0) * log_stdv - sd_lambda * b * b / ev[i].stdv) / static_cast< float >(2.0);
		
		m[i*n_states+j].alpha += log_normal_pdf+log_invgauss_pdf;

            }
        }
	
}


#endif
