#ifndef __DATA_STRUCTURE_HPP
#define __DATA_STRUCTURE_HPP


#include<iostream>
struct Matrix_Entry{
	float alpha;
	unsigned beta;
};
struct event{
	float mean, stdv,start, length;
	
};
struct pore_model_state{
	float level_mean, level_stdv,sd_mean, sd_stdv;
};


#endif
