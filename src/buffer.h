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
#include <stdio.h>
#include <omp.h>

// para:
//		filename
//		nrow/ncol
//      output
template <class T> 
void read_in_buffer(std::string filename, int nrow, int ncol, std::vector<T> ouput){
	/*
	READING TO BUFFER;
	*/
	long size_len = (long)nrow*((long)ncol); // verify size_len == size
	// buffer read, type = float
	size_t mul = sizeof(float) / sizeof(unsigned char);
	FILE *fp;
	fp = fopen(filename.c_str(), "rb");// open as binary
	fseek(fp, 0, SEEK_END);
	long size = ftell(fp);         /*calc the size needed*/
	fseek(fp, 0, SEEK_SET);
	unsigned char * buffer = (unsigned *)malloc(size); /*allocalte space on heap*/
	float * buffer_float = (float *)buffer;// recasting to float
	if (fp == NULL){ /*ERROR detection if file == empty*/
		printf("Error: There was an Error reading the file %s \n", path);
		exit(1);
	}
	else if (fread(&buffer, sizeof(*buffer), size, fp) != size){ /* if count of read bytes != calculated size of .bin file -> ERROR*/
		printf("Error: There was an Error reading the file %s - %d\n", path, r);
		exit(1);
	}
	else{
		int i;
		for (i = 0; i<size/mul; i++){
			std::cout << buffer_float[i];
		}
	}
	free(buffer);
	fclose(fp);
}
#endif