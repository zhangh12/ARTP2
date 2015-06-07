
#include <R.h>
#include <stdio.h>
#include <iostream>

using namespace std;

typedef unsigned char uint8;

extern "C" {

void ReadBED(char **input_bed_file, int *input_nsub, int *input_nsnp, 
int *input_nsel, int *sel_snp_id, int *geno){
  
  int nsub = *input_nsub;
  //int nsnp = *input_nsnp;
  int nsel = *input_nsel;
  
  FILE *bed_file = fopen(*input_bed_file, "rb");
  if(!bed_file){
    error("Could not open %s", *input_bed_file);
  }
  
  unsigned char magic[3];
  if(fread(magic, 1, 3, bed_file) != 3){
    error("Failed to read first three bytes");
  }
  
  if(magic[0] != 0x6c || magic[1] != 0x1b){
    error("%s is not a binary PED file", *input_bed_file);
  }
  
  if(magic[2] != 0x01){
    error("%s is not in SNP-major mode", *input_bed_file);
  }
  
  uint8 PROBE1[4] = {0x01, 0x04, 0x10, 0x40};
  uint8 PROBE2[4] = {0x02, 0x08, 0x20, 0x80};
  
  int nblock = (int) ceil((double) nsub / 4);
  
  int offset = -1;
  for(int i = 0; i < nsel; ++i){
    int ns = 0;
    if(i == 0){
      ns = sel_snp_id[i] - 1;
    }else{
      ns = sel_snp_id[i] - sel_snp_id[i - 1] - 1;
    }
    for(int j = 0; j < ns; ++j){
      fseek(bed_file, nblock, SEEK_CUR);
    }
    
    int sid = 0;
    for(int k = 0; k < nblock; ++k){
      uint8 b = fgetc(bed_file);
      if(feof(bed_file)){
        break;
      }
      for(int l = 0; l < 4; ++l){
        uint8 h1 = b & PROBE1[l];
        uint8 h2 = b & PROBE2[l];
        
        ++offset;
        if(h1 && !h2){// missing
          geno[offset] = -1;
        }else if(!h1 && !h2){// 0/0
          geno[offset] = 0;
        }else if(!h1 && h2){// 0/1
          geno[offset] = 1;
        }else{// 1/1
          geno[offset] = 2;
        }
        
        ++sid;
        if(sid >= nsub){
          break;
        }
      }
    }
  }
  
  fclose(bed_file);
  
}

}
