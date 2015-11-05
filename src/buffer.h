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
	size_t mul = sizeof(float) / sizeof(char); // mul >= 2
	const long size_len = (long)nrow*((long)ncol)*mul; // verify size_len == size
	// buffer read, type = float
    char * buffer = (char *)malloc(sizeof(char)*size_len);
	std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);
	fin.read(buffer, size_len);
	float * buffer_float = (float *)buffer;  //TODO: recasting char* to float *
	for (int i = 0; i < nrow; i++){
		std::cout << "\n[" << i << ",]:\t";
		for (int j = 0; j < ncol; j++)
			std::cout  << buffer_float[i*ncol + j]<<"\t";
	}
	std::cout << "\n";
	fin.close();
	free(buffer);
}
#endif

/*TEST FILE LIST
>>>>>>>>>>>>>  testbuffer1.cpp >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <omp.h>
#include <vector>

#include <time.h>

#include <sys/time.h>
#include "../ARTP3/src/buffer.h"
class StopWatch {
	timeval     started;
	std::string msg;

public:
	StopWatch(const std::string& m) : msg(m)
	{
		gettimeofday(&started, NULL);
	}

	~StopWatch() {
		std::cerr << msg << " " << usecs() << " µsecs" << std::endl;
	}

	unsigned usecs() const {
		timeval t2;
		gettimeofday(&t2, NULL);

		return (t2.tv_sec - started.tv_sec) * 1000000 + (t2.tv_usec - started.tv_usec);
	}
};
using namespace std;

int main()
{
	int nrow = 10;
	int ncol = 2;
	std::string filename("raw.dta");
	{
		StopWatch w1("Write file in");
		ofstream f1("raw.dta", ios::binary | ios::out);
		float f = 3.14159f;
		if (f1.is_open()) {
			for (int i = 0; i < nrow; i++) {
				for (int j = 0; j < ncol; j++) {
					float ff = f*(i + j);
					f1.write(reinterpret_cast<const char*>(&ff), sizeof(f));
				}
			}
			f1.close();
		}
	}
	{// TEST read buffer
	std:vector<float> U;
		read_in_buffer(filename, 10, 2, U);

	}
	return 0;
}
>>>>>>>>>>>>>EOF  testbuffer1.cpp >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/