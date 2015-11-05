#ifndef BUFFER_BIN_ARTPS
#define BUFFER_BIN_ARTPS
/*
Task:
Buffer reading binary file,
given file size (nrow ncol)
all data are float

.815
artp3.cpp:
- artp3_chr    L672-683
- adajoint_chr L924-936
- artp3        L1059-1069

*/

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// para:
//		filename
//		nrow/ncol
//      output
template <class T> 
void read_in_buffer(std::string filename, int nrow, int ncol, std::vector<T> ouput){
	/*
	TODO: READING TO BUFFER;
	Verified: testbuffer1.cpp
	*/
	size_t mul = sizeof(float) / sizeof(unsigned char);
	const long size_len = (long)nrow*((long)ncol)*mul; // verify size_len == size
	// buffer read, type = float
    char * buffer = (char *)malloc(sizeof(char)*size_len);
	std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);
	fin.read(buffer, size_len);
	float * buffer_float = (float *)buffer;
	for (int i = 0; i < nrow; i++){
		std::cout << "\n[" << i << ",]:\t";
		for (int j = 0; j < ncol; j++)
			std::cout  << buffer_float[i*ncol + j]<<"\t";
	}
	std::cout << "\n";

}
#endif