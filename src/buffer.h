#ifndef BUFFER_BIN_ARTPS
#define BUFFER_BIN_ARTPS

//#define PRINT_DEB
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
int segment_read(char *buff, const int len, const int count) {
	return 1;
}

template <class T> 
void read_in_buffer(std::string filename, int nrow, int ncol, std::vector<T>& output) {

	/*
	TODO: READING TO BUFFER;
	Verified: testbuffer1.cpp
	*/
	std::cout << "Read into buffer" << std::endl;
	size_t mul = sizeof(float) / sizeof(char); // mul >= 2
	const long size_len = (long)nrow*((long)ncol)*mul; // verify size_len == size
	size_t sizeA = nrow*ncol;
	// buffer read, type = float
	char * buffer = (char *)malloc(sizeof(char)*size_len);
	if (buffer == NULL) {
		std::cout << "Out of Memory." << std::endl;
		exit(1);
	}
	std::ifstream fin(filename.c_str(), std::ios::in | std::ios::binary);
	fin.read(buffer, size_len);
	float * buffer_float = (float *)buffer;  //TODO: recasting char* to float *
#ifdef PRINT_DEB
	for (int i = 0; i < nrow; i++) {
		std::cout << "\n[" << i << ",]:\t";
		for (int j = 0; j < ncol; j++)
			std::cout << buffer_float[i*ncol + j] << "\t";
	}
	std::cout << "\n";
#endif
	fin.close();
	std::cout << "Read into OMP" << std::endl;
		int nthreads;
		int i;
#pragma omp parallel
		{
			nthreads = omp_get_num_threads();
			std::cout << "Thread:" << nthreads <<std::endl;
#pragma omp parallel for num_threads(nthreads) private(i)  shared(output)
			for (i = 1; i < sizeA; i++) {
				output[i] = buffer_float[i];
				
			}
		}
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
int nrow = 12345679;
int ncol = 7;
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
StopWatch ww2("Read file in OMP");
std:vector<float> U;
U.resize(nrow*ncol);
read_in_buffer(filename, nrow, ncol, U);
//for (int i = 0; i < nrow*ncol; i++) cout << U[i] << '\t';
}
return 0;
}
>>>>>>>>>>>>>EOF  testbuffer1.cpp >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/