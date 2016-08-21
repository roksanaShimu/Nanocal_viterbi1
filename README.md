# Nanocal_viterbi1


new.cu //this is the main file

Data_Structure.hpp // this contains the data stuctures (struct, class  etc) used in the program

Predefine_Values.hpp // lists all predefine values like block size, number of threads etc.

Initialization.hpp //In the main program, there are some arrays / vectors containing neighbors' state_number , probabilities, mean,... etc. 			//These values should come from nanocall. But as I am working on a fixed block (a part of "fill" function from Viterbi.hpp), I 			//have generated random numbers for the arrays'.  Initialization.hpp contains the code for generating the random numbers. In 			//future, when I add the parallel code to the main nanocall program, I won't need these codes.

Serial_Implementation.hpp //serial implementation of the block I am working on. copied from the nanocall program

Parallel_Implementation_1i.hpp // parallel implementation of approach 1 (picture is attached in the work update file dropbox)

Device_Memory_Allocation.hpp // contains functions used to allocate different types of memory.  

Parallel_Implementation_1ii.hpp // parallel implementation of approach 1 (picture is attached in the work update file dropbox)
			 // the difference between Parallel_Implementation_1i.hpp and Parallel_Implementation_1ii.hpp is the memory hierarchy
