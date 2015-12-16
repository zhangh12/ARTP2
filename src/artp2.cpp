#include "util.h"


struct STAT{
  
  float stat;
  int id;
  
  bool operator < (const STAT& st) const{ // ascending
    return (stat < st.stat);
  }
  
  STAT(float st, int i){
    stat = st;
    id = i;
  }
  
};

bool descending(const STAT &a, const STAT &b){ return a.stat > b.stat; }
bool ascending(const STAT &a, const STAT &b){ return a.stat < b.stat; }

STAT STAT0(.0f, -1);

struct MINP{
  
  int stat;
  int id;
  
  bool operator < (const MINP &mp) const{
    return (stat < mp.stat);
  }
  
  MINP(int st, int i){
    stat = st;
    id = i;
  }
  
};

MINP MINP0(0, -1);

typedef vector<STAT> VecStat;
typedef vector<MINP> VecMinp;
typedef vector<float> fvec;
typedef vector<vector<float> > fmat;
typedef vector<double> dvec;
typedef vector<vector<double> > dmat;
typedef vector<int> ivec;
typedef vector<vector<int> > imat;
typedef vector<string> svec;
typedef vector<vector<string> > smat;

void load_U(double *R_vU, fmat &U, int nsnp){
  
  U = fmat (nsnp, fvec (nsnp, .0f));
  int k = -1;
  for(int i = 0; i < nsnp; ++i){
    for(int j = 0; j < nsnp; ++j){
      ++k;
      U[i][j] = R_vU[k];
    }
  }
  
}

void load_score0(double *R_score0, fvec &score0, int nsnp){
  
  score0 = fvec (nsnp, .0f);
  for(int i = 0; i < nsnp; ++i){
    score0[i] = R_score0[i];
  }
  
}

void load_cov(double *R_vV, fmat &V, int nsnp){
  
  V = fmat (nsnp, fvec (nsnp, .0f));
  int k = -1;
  for(int i = 0; i < nsnp; ++i){
    for(int j = 0; j < nsnp; ++j){
      ++k;
      V[i][j] = R_vV[k];
    }
  }
  
}

void load_sigma2(const fmat &V, fvec &sigma2){
  
  int nsnp = V.size();
  sigma2 = fvec (nsnp, .0f);
  for(int i = 0; i < nsnp; ++i){
    sigma2[i] = V[i][i];
  }
}

void load_gene_idx(int *R_vgene_idx, int *R_gene_start, int *R_gene_end, 
imat &gene_idx, int ngene){
  
  gene_idx = imat (ngene);
  for(int i = 0; i < ngene; ++i){
  	ivec idx;
  	for(int j = R_gene_start[i]; j <= R_gene_end[i]; ++j){
  		idx.push_back(R_vgene_idx[j - 1] - 1);
  	}
  	gene_idx[i] = idx;
  }
  
}

void load_gene_cutpoint(int *R_vgene_cutpoint, int *R_gene_cutpoint_start, int *R_gene_cutpoint_end, 
imat &cutpoint, int ngene){
  
  cutpoint = imat (ngene);
  for(int i = 0; i < ngene; ++i){
  	ivec cp;
  	for(int j = R_gene_cutpoint_start[i]; j <= R_gene_cutpoint_end[i]; ++j){
  		cp.push_back(R_vgene_cutpoint[j - 1] - 1);
  	}
    cutpoint[i] = cp;
  }
  
}


//gene_idx is for a single gene (otherwise it should be imat rather than ivec)
void extract_score(fvec &S, const fvec &score0, 
const ivec &gene_idx){
  
	int ns = gene_idx.size();//number of SNPs in a gene
		
	S = fvec(ns);
	for(int j = 0; j < ns; ++j){
		int J = gene_idx[j];
		S[j] = score0[J];
	}
	
}

void extract_cov(fmat &Sigma, const fmat &V, 
const ivec &gene_idx){
	
	int ns = gene_idx.size();//number of SNPs in a gene
		
	Sigma = fmat(ns, fvec(ns));
	for(int j = 0; j < ns; ++j){
		int J = gene_idx[j];
		Sigma[j][j] = V[J][J];
		for(int k = j + 1; k < ns; ++k){
			Sigma[j][k] = V[J][gene_idx[k]];
			Sigma[k][j] = Sigma[j][k];
		}
	}
	
}


void search1(fvec &stat, ivec &sel_id, int &marg_id, const fvec &S, const fmat &Sigma, const int &mc){
  
  int ns = S.size();
  
	//find best marginal SNP
	int id = -1;
	double max_stat = -1.0;
	for(int j = 0; j < ns; ++j){
		double s = S[j] * S[j] / Sigma[j][j];
		if(s > max_stat){
			max_stat = s;
			id = j;
		}
	}
	
  stat = fvec(1, max_stat);
	
	assert(id >= 0 && id < ns);
	sel_id = ivec(1, id);
  marg_id = id;
  
	if(ns == 1){
		return;
	}
	
	dmat inv_Sigma22(mc + 1, dvec(mc + 1, .0));
	inv_Sigma22[0][0] = 1.0 / Sigma[id][id];
	
	dvec S2(mc + 1, .0);
	S2[0] = S[id];
	
	vector<bool> used(ns, false);
	used[id] = true;
	
	double w = max_stat;
	
	for(int j = 1; j <= mc; ++j){//j is the number of SNPs being selected
		
		dvec v(j, .0); // v = (Sigma22)^-1 S2
		for(int i = 0; i < j; ++i){
			v[i] = .0;
			for(int k = 0; k < j; ++k){
				v[i] += inv_Sigma22[i][k] * S2[k];
			}
		}
		
		max_stat = -1.0;
		double max_qf2 = -1.0;
		id = -1;
		for(int k = 0; k < ns; ++k){
			if(used[k]){
				continue;
			}
			
			double S1 = S[k];
			double qf1 = .0; //qf1 = sigma12 (Sigma22)^-1 S2
			double qf2 = .0; //qf2 = sigma12 (Sigma22)^-1 sigma21
			for(int l = 0; l < j; ++l){
				qf1 += Sigma[k][sel_id[l]] * v[l];
				double tmp = .0; // tmp = the l-th column of sigma12 (Sigma22)^-1
				for(int t = 0; t < j; ++t){
					tmp += Sigma[k][sel_id[t]] * inv_Sigma22[t][l];
				}
				qf2 += tmp * Sigma[sel_id[l]][k];
			}
			
			
			double stat = w + (S1 * S1 - 2.0 * S1 * qf1 + qf1 * qf1)/(Sigma[k][k] - qf2);
			if(stat > max_stat){
				max_stat = stat;
				id = k;
				max_qf2 = qf2;
			}
			
		}
		
		w = max_stat;
		stat.push_back(max_stat);
    
		sel_id.push_back(id);
		used[id] = true;
		S2[j] = S[id];
		inv_Sigma22[j][j] = 1.0 / (Sigma[id][id] - max_qf2);
		dvec tmp(j, .0);
		for(int l = 0; l < j; ++l){
			inv_Sigma22[j][l] = .0;
			for(int t = 0; t < j; ++t){
				inv_Sigma22[j][l] += Sigma[id][sel_id[t]] * inv_Sigma22[t][l];
			}
			inv_Sigma22[j][l] *= -inv_Sigma22[j][j];
			inv_Sigma22[l][j] = inv_Sigma22[j][l];
		}
		
		for(int l = 0; l < j; ++l){
			for(int t = 0; t < j; ++t){
				inv_Sigma22[l][t] +=  inv_Sigma22[l][j] * inv_Sigma22[t][j] / inv_Sigma22[j][j];
			}
		}
		
	}
	
}

void search2(fvec &stat, ivec &sel_id, int &marg_id, const fvec &S, const fmat &Sigma, const int &mc){
  
  int ns = S.size();
  
	//find best marginal SNP
	int id = -1;
	double max_stat = -1.0;
	for(int j = 0; j < ns; ++j){
		double s = S[j] * S[j] / Sigma[j][j];
		if(s > max_stat){
			max_stat = s;
			id = j;
		}
	}
	
  stat = fvec(1, max_stat);
  assert(id >= 0 && id < ns);
  sel_id = ivec(1, id);
  marg_id = id;
  
	if(ns == 1){
		return;
	}
  
  //find best pair
  int id1 = -1;
  int id2 = -1;
  max_stat = -1.0;
  for(int j = 0; j < ns - 1; ++j){
    for(int k = j + 1; k < ns; ++k){
      double a = Sigma[j][j];
      double b = Sigma[j][k];
      double d = Sigma[k][k];
      double stat = (d * S[j] * S[j] - 2 * b * S[j] * S[k] + a * S[k] * S[k]) / (a * d - b * b);
      if(stat > max_stat){
        max_stat = stat;
        id1 = j;
        id2 = k;
      }
    }
  }
  
  stat.push_back(max_stat);
  assert(id1 >= 0 && id1 < ns);
  assert(id2 >= 0 && id2 < ns);
  sel_id = ivec(1, id1);
  sel_id.push_back(id2);
  
  if(ns == 2){
    return;
  }
	
	dmat inv_Sigma22(mc + 1, dvec(mc + 1, .0));
  double a = Sigma[id1][id1];
  double b = Sigma[id1][id2];
  double d = Sigma[id2][id2];
  double e = a * d - b * b;
	inv_Sigma22[0][0] = d / e;
  inv_Sigma22[0][1] = -b / e;
  inv_Sigma22[1][0] = inv_Sigma22[0][1];
  inv_Sigma22[1][1] = a / e;
	
	dvec S2(mc + 1, .0);
	S2[0] = S[id1];
  S2[1] = S[id2];
	
	vector<bool> used(ns, false);
	used[id1] = true;
  used[id2] = true;
	
	double w = max_stat;
	
	for(int j = 2; j <= mc; ++j){//j is the number of SNPs being selected
		
		dvec v(j, .0); // v = (Sigma22)^-1 S2
		for(int i = 0; i < j; ++i){
			v[i] = .0;
			for(int k = 0; k < j; ++k){
				v[i] += inv_Sigma22[i][k] * S2[k];
			}
		}
		
		max_stat = -1.0;
		double max_qf2 = -1.0;
		id = -1;
		for(int k = 0; k < ns; ++k){
			if(used[k]){
				continue;
			}
			
			double S1 = S[k];
			double qf1 = .0; //qf1 = sigma12 (Sigma22)^-1 S2
			double qf2 = .0; //qf2 = sigma12 (Sigma22)^-1 sigma21
			for(int l = 0; l < j; ++l){
				qf1 += Sigma[k][sel_id[l]] * v[l];
				double tmp = .0; // tmp = the l-th column of sigma12 (Sigma22)^-1
				for(int t = 0; t < j; ++t){
					tmp += Sigma[k][sel_id[t]] * inv_Sigma22[t][l];
				}
				qf2 += tmp * Sigma[sel_id[l]][k];
			}
			
			
			double stat = w + (S1 * S1 - 2.0 * S1 * qf1 + qf1 * qf1)/(Sigma[k][k] - qf2);
			if(stat > max_stat){
				max_stat = stat;
				id = k;
				max_qf2 = qf2;
			}
			
		}
		
		w = max_stat;
		stat.push_back(max_stat);
    
		sel_id.push_back(id);
		used[id] = true;
		S2[j] = S[id];
		inv_Sigma22[j][j] = 1.0 / (Sigma[id][id] - max_qf2);
		dvec tmp(j, .0);
		for(int l = 0; l < j; ++l){
			inv_Sigma22[j][l] = .0;
			for(int t = 0; t < j; ++t){
				inv_Sigma22[j][l] += Sigma[id][sel_id[t]] * inv_Sigma22[t][l];
			}
			inv_Sigma22[j][l] *= -inv_Sigma22[j][j];
			inv_Sigma22[l][j] = inv_Sigma22[j][l];
		}
		
		for(int l = 0; l < j; ++l){
			for(int t = 0; t < j; ++t){
				inv_Sigma22[l][t] +=  inv_Sigma22[l][j] * inv_Sigma22[t][j] / inv_Sigma22[j][j];
			}
		}
		
	}
	
}

#if __PARALLEL__
void rnorm(drand48_data &buf, int n, fvec &v){
	
	v = fvec (n);
	
	int m = -1;
	if(n % 2){
		m = n/2 + 1;
	}else{
		m = n/2;
	}
	
  int k = -1;
	for(int i = 0; i < m; ++i){
		float u1, u2, s;
		double r1, r2;
		do{
			drand48_r(&buf, &r1);
			drand48_r(&buf, &r2);
			u1 = 2.0f * (float) r1 - 1.0f;
			u2 = 2.0f * (float) r2 - 1.0f;
			s = u1 * u1 + u2 * u2;
		}while(s >= 1.0f || (s < 1e-12 && s > -1e-12));
		
		s = (float) sqrt((-2.0 * log(s)) / s);
    ++k;
		v[k] = u1 * s;
    ++k;
    if(k < n){
      v[k] = u2 * s;
    }
	}
	
}
#else
void rnorm(int n, fvec &v){
  
	v = fvec (n);
	
	int m = -1;
	if(n % 2){
		m = n/2 + 1;
	}else{
		m = n/2;
	}
	
  int k = -1;
	for(int i = 0; i < m; ++i){
		float u1, u2, s;
		double r1, r2;
		do{
			r1 = rand() / ((double) RAND_MAX);
			r2 = rand() / ((double) RAND_MAX);
			u1 = 2.0f * (float) r1 - 1.0f;
			u2 = 2.0f * (float) r2 - 1.0f;
			s = u1 * u1 + u2 * u2;
		}while(s >= 1.0f || (s < 1e-12 && s > -1e-12));
		
		s = (float) sqrt((-2.0 * log(s)) / s);
    ++k;
		v[k] = u1 * s;
    ++k;
    if(k < n){
      v[k] = u2 * s;
    }
	}
	
}
#endif

void load_group_id(int *R_group_id, ivec &group_id, int ngene){
  
  group_id = ivec(ngene);
  for(int i = 0; i < ngene; ++i){
    group_id[i] = R_group_id[i] - 1;
  }
  
}
  
void load_gene_id(int *R_gene_id, ivec &gene_id, int ngene){
  
  gene_id = ivec(ngene);
  for(int i = 0; i < ngene; ++i){
    gene_id[i] = R_gene_id[i] - 1;
  }
  
}

void load_pathway_cutpoint(int *R_pathway_cutpoint, ivec &pathway_cutpoint, int ncp){
  
  pathway_cutpoint = ivec(ncp);
  for(int i = 0; i < ncp; ++i){
    pathway_cutpoint[i] = R_pathway_cutpoint[i] - 1;
  }
  
}

void read_in_buffer(const string &file, const int &nperm, const int &ncp, const int &nthread, vector<VecStat> &stat){
  
  stat = vector<VecStat>(ncp, VecStat (nperm + 1, STAT0));
  
  FILE * gFile = fopen(file.c_str(), "rb");
  if(gFile == NULL){
    error("Cannot open gene output file");
  }
  
  //check file size
  fseek(gFile, 0, SEEK_END);
  long lsize = ftell(gFile);
  rewind(gFile);
  assert(lsize == (long) (nperm + 1L) * (long) (ncp) * sizeof(float)/sizeof(char));
  
  char *buffer = new char[lsize];
  if(buffer == NULL){
    error("Out of memory in loading gene output file");
  }
  
  size_t status = fread(buffer, 1, lsize, gFile);
  if(status != lsize){
    error("Gene output file might be modified by other jobs in queue. Please check options$id.str and options$out.dir");
  }
  
  fclose(gFile);
  
  float *buffer_float = (float *) buffer;
  #if __PARALLEL__
  #pragma omp parallel num_threads(nthread)
  {
    #pragma omp for
    #endif
    for(int i = 0; i < nperm + 1; ++i){
      for(int j = 0; j < ncp; ++j){
        stat[j][i].stat = buffer_float[i * ncp + j];
        stat[j][i].id = i;
      }
    }
  #if __PARALLEL__
  }
  #endif
  
  buffer_float = NULL;
  delete[] buffer;
  
}

extern "C" {

void artp2_chr(char **R_file_prefix, int *R_method,
int *R_nperm, int *R_seed, 
int *R_nthread, int *R_nsnp, int *R_ngene, 
double *R_vU, double *R_score0, double *R_vV, 
int *R_vgene_idx, int *R_gene_start, int *R_gene_end, 
int *R_vgene_cutpoint, 
int *R_gene_cutpoint_start, int *R_gene_cutpoint_end, 
double *R_gene_pval, int *R_arr_rank, 
int *R_sel_id, int *R_marg_id){
  
  int len_file_prefix = strlen(*R_file_prefix);
  char *file_prefix = new char[len_file_prefix + 1];
  file_prefix[0] = '\0';
  strcat(file_prefix, *R_file_prefix);
  
  int method = *R_method;
  assert(method == 3);
  if(method == 3){
    ;
  }
  
  int nperm = *R_nperm;
  int seed = *R_seed;
  int nthread = *R_nthread;
  int nsnp = *R_nsnp;
  int ngene = *R_ngene;
  
  #if __NO_PARALLEL__
  assert(nthread == 1);
	#endif
	
  fvec score0;
  fvec sigma2;
  fmat U;
  fmat V;
  load_score0(R_score0, score0, nsnp);
  load_cov(R_vV, V, nsnp);
  load_sigma2(V, sigma2);
  load_U(R_vU, U, nsnp);
  
  imat gene_idx; // index of SNPs in a gene
  
  load_gene_idx(R_vgene_idx, R_gene_start, R_gene_end, 
  gene_idx, ngene);
  
  imat cutpoint;
  load_gene_cutpoint(R_vgene_cutpoint, R_gene_cutpoint_start, R_gene_cutpoint_end, 
  cutpoint, ngene);
  
  string fprefix (file_prefix);
  svec gene_out (ngene, fprefix);
  for(int g = 0; g < ngene; ++g){
    ostringstream gid;
    gid << g;
    gene_out[g] = gene_out[g] + string("GID.") + gid.str() + string(".bin");
  }
  
  // write obs statistics for all genes
  imat sel_id(ngene);
  ivec marg_id(ngene);
  for(int g = 0; g < ngene; ++g){
  	fstream gout(gene_out[g].c_str(), ios::out | ios::binary);
  	if(!gout){
  		error("Fail to write observed statistics to file");
  	}
  	int ns = gene_idx[g].size();
    int ncp = cutpoint[g].size();
    int max_cutpoint = cutpoint[g][ncp - 1];
    fvec s (ns, .0f);
    VecStat vs (ns, STAT0);
  	for(int j = 0; j < ns; ++j){
  		s[j] = score0[gene_idx[g][j]];
      s[j] = pchisq(s[j] * s[j] / sigma2[gene_idx[g][j]], 1, false, true);
      vs[j].stat = -s[j];
      vs[j].id = j;
  	}
    
    sort(vs.begin(), vs.end(), descending);
    marg_id[g] = vs[0].id;
    for(int j = 0; j < ns; ++j){
      sel_id[g].push_back(vs[j].id);
    }
    
    sort(s.begin(), s.end());
    for(int j = 1; j <= max_cutpoint; ++j){
      s[j] += s[j - 1];
    }
    
    for(int k = 0; k < ncp; ++k){
      float u = -s[cutpoint[g][k]];
      gout.write((char*)(&u), sizeof(u));
    }
  	gout.close();
  }
  
  int i_sel_id = -1;
  for(int g = 0; g < ngene; ++g){
    R_marg_id[g] = gene_idx[g][marg_id[g]] + 1;
    for(int k = 0; k < (int) sel_id[g].size(); ++k){
      ++i_sel_id;
      R_sel_id[i_sel_id] = gene_idx[g][sel_id[g][k]] + 1;
    }
    int nn = gene_idx[g].size() - sel_id[g].size();
    while(nn){
      ++i_sel_id;
      R_sel_id[i_sel_id] = -1;
      --nn;
    }
  }
  
  int ngap = min(10000, nperm);
  int nblock = nperm / ngap;
  
	#if __NO_PARALLEL__
  srand(seed);
	#endif
	
  for(int b = 0; b < nblock; ++b){
  	fmat null(ngap, fvec (nsnp, .0f));
		#if __PARALLEL__
  	drand48_data buf;
		#endif
  	// compute null statistics
  	#if __PARALLEL__
  	#pragma omp parallel num_threads(nthread) private(buf)
  	{
  		srand48_r(seed + b * nthread + omp_get_thread_num(), &buf);
  		#pragma omp for
  	#endif
  		for(int i = 0; i < ngap; ++i){
	  		fvec rn;
				#if __PARALLEL__
	  		rnorm(buf, nsnp, rn);
				#else
				rnorm(nsnp, rn);
				#endif
	  		for(int j = 0; j < nsnp; ++j){
	  			null[i][j] = .0f;
	  			for(int k = 0; k < nsnp; ++k){
	  				null[i][j] += rn[k] * U[k][j];
	  			}
	  			null[i][j] = null[i][j] * null[i][j] / sigma2[j];
	  			null[i][j] = pchisq(null[i][j], 1, false, true);
	  		}
	  	}
		#if __PARALLEL__
	  }
		#endif
	  
	  // write null statistics to local files (per gene)
		#if __PARALLEL__
	  #pragma omp parallel num_threads(min(nthread, ngene))
	  {
	  	#pragma omp for
	  #endif
	  	for(int g = 0; g < ngene; ++g){
	  		ofstream gout;
	  		gout.open(gene_out[g].c_str(), ios::out | ios::binary | ios::app);
	  		if(!gout){
	  			error("Fail to write null statistics to file");
	  		}
        int ns = gene_idx[g].size();
        int ncp = cutpoint[g].size();
        int max_cutpoint = cutpoint[g][ncp - 1];
	  		for(int i = 0; i < ngap; ++i){
	  			fvec s(ns, .0f);
	  			for(int j = 0; j < ns; ++j){
	  				s[j] = null[i][gene_idx[g][j]];
	  			}
          
          sort(s.begin(), s.end());
          for(int j = 1; j <= max_cutpoint; ++j){
            s[j] += s[j - 1];
          }
          
          for(int k = 0; k < ncp; ++k){
            float u = -s[cutpoint[g][k]];
            gout.write((char*)(&u), sizeof(u));
          }
	  		}
	  		gout.close();
	  	}
		#if __PARALLEL__
	  }
		#endif
	  fmat().swap(null);
  }
  
  // read null statistics (per gene)
  int irk = -1;
  for(int g = 0; g < ngene; ++g){
  	int ncp = cutpoint[g].size();
    
    vector<VecStat> stat(ncp, VecStat (nperm + 1, STAT0));
    fstream gin(gene_out[g].c_str(), ios::in | ios::binary);
    for(int i = 0; i < nperm + 1; ++i){
  		for(int j = 0; j < ncp; ++j){
  			float s = .0f;
  			gin.read((char*)(&s), sizeof(s));
  			stat[j][i].stat = s;
  			stat[j][i].id = i;
  		}
  	}
    gin.close();
    /*
  	vector<VecStat> stat;
    read_in_buffer(gene_out[g], nperm, ncp, nthread, stat);
    */
    
    if(remove(gene_out[g].c_str())){
      error("Cannot delete gene output file");
    }
  	
  	imat arr_rank(ncp, ivec (nperm + 1, 0));
		#if __PARALLEL__
  	#pragma omp parallel num_threads(min(ncp, nthread))
  	{
      #pragma omp for
    #endif
    	for(int j = 0; j < ncp; ++j){
  			sort(stat[j].begin(), stat[j].end(), descending);
  			for(int i = 0; i < nperm + 1; ++i){
  				int id = stat[j][i].id;
  				arr_rank[j][id] = i;
  			}
  			VecStat().swap(stat[j]);
  		}
		#if __PARALLEL__
  	}
		#endif
  	
  	vector<VecStat>().swap(stat);
    
    VecMinp gene_min_p (nperm + 1, MINP0);
    ivec subsum(nthread, 0);
    ivec subtie(nthread, 0);
    int m = nperm + 1;
    for(int j = 0; j < ncp; ++j){
      ++irk;
      R_arr_rank[irk] = arr_rank[j][0];
      if(arr_rank[j][0] < m){
        m = arr_rank[j][0];
      }
    }
    gene_min_p[0].stat = m;
    gene_min_p[0].id = 0;
    
    #if __PARALLEL__
    #pragma omp parallel num_threads(nthread)
    {
      #pragma omp for
    #endif
      for(int i = 1; i < nperm + 1; ++i){
      	#if __PARALLEL__
        int tid = omp_get_thread_num();
        #else
        int tid = 0;
        #endif
        int m = nperm + 1;
        for(int j = 0; j < ncp; ++j){
          if(arr_rank[j][i] < m){
            m = arr_rank[j][i];
          }
        }
        gene_min_p[i].stat = m;
        gene_min_p[i].id = i;
        if(gene_min_p[i].stat < gene_min_p[0].stat){
          subsum[tid] += 1;
        }else if(gene_min_p[i].stat == gene_min_p[0].stat){
          subtie[tid] += 1;
        }else{
          ;
        }
      }
    #if __PARALLEL__
    }
    #endif
    
    R_gene_pval[g] = 1.0;
    int rep = 0;
    for(int t = 0; t < nthread; ++t){
      R_gene_pval[g] += subsum[t];
      rep += subtie[t];
    }
    R_gene_pval[g] += rep / 2.0;
    R_gene_pval[g] /= nperm + 1;
    
    sort(gene_min_p.begin(), gene_min_p.end());
    
    int *gene_p_stat = new int[nperm + 1];
    if(gene_p_stat == NULL){
    	error("Out of memory in artp2_chr");
    }
    
    #if __PARALLEL__
    #pragma omp parallel num_threads(nthread)
    {
    	#pragma omp for
    #endif
	    for(int i = 0; i < nperm + 1; ++i){
	    	int id = gene_min_p[i].id;
	    	gene_p_stat[id] = i;
	    }
	  #if __PARALLEL__
	  }
	  #endif
    VecMinp().swap(gene_min_p);
    
    fstream gout(gene_out[g].c_str(), ios::out | ios::binary);
    if(!gout){
      error("Fail to write gene statistics to file");
    }
    
    gout.write((char*)(gene_p_stat), sizeof(int) * (nperm + 1));
    delete[] gene_p_stat;
    gout.close();
    
  }
  
  delete[] file_prefix;
  
}

void adajoint_chr(char **R_file_prefix, int *R_method,
int *R_nperm, int *R_seed, 
int *R_nthread, int *R_nsnp, int *R_ngene, 
double *R_vU, double *R_score0, double *R_vV, 
int *R_vgene_idx, int *R_gene_start, int *R_gene_end, 
int *R_vgene_cutpoint, 
int *R_gene_cutpoint_start, int *R_gene_cutpoint_end, 
double *R_gene_pval, int *R_arr_rank, 
int *R_sel_id, int *R_marg_id){
  
  int len_file_prefix = strlen(*R_file_prefix);
  char *file_prefix = new char[len_file_prefix + 1];
  file_prefix[0] = '\0';
  strcat(file_prefix, *R_file_prefix);
  
  int method = *R_method;
  assert(method == 1 || method == 2);
  
  int nperm = *R_nperm;
  int seed = *R_seed;
  int nthread = *R_nthread;
  int nsnp = *R_nsnp;
  int ngene = *R_ngene;
  
  #if __NO_PARALLEL__
  assert(nthread == 1);
  #endif
  
  fvec score0;
  fmat V;
  fmat U;
  load_score0(R_score0, score0, nsnp);
  load_cov(R_vV, V, nsnp);
  load_U(R_vU, U, nsnp);
  
  imat gene_idx; // index of SNPs in a gene
  
  load_gene_idx(R_vgene_idx, R_gene_start, R_gene_end, 
  gene_idx, ngene);
  
  imat cutpoint;
  load_gene_cutpoint(R_vgene_cutpoint, R_gene_cutpoint_start, R_gene_cutpoint_end, 
  cutpoint, ngene);
  
  string fprefix (file_prefix);
  svec gene_out (ngene, fprefix);
  for(int g = 0; g < ngene; ++g){
    ostringstream gid;
    gid << g;
    gene_out[g] = gene_out[g] + string("GID.") + gid.str() + string(".bin");
  }
  
  // write obs statistics for all genes
  imat sel_id(ngene);
  ivec marg_id(ngene);
  for(int g = 0; g < ngene; ++g){
    fstream gout(gene_out[g].c_str(), ios::out | ios::binary);
    if(!gout){
      error("Fail to write observed statistics to file");
  	}
  	
  	fvec S;
  	fmat Sigma;
  	extract_score(S, score0, gene_idx[g]);
  	extract_cov(Sigma, V, gene_idx[g]);
  	fvec s;
    int ncp = cutpoint[g].size();
    int mc = cutpoint[g][ncp - 1];
    
    if(method == 1){
      search1(s, sel_id[g], marg_id[g], S, Sigma, mc);
    }else{//assert(method == 2)
      search2(s, sel_id[g], marg_id[g], S, Sigma, mc);
    }
    
    for(int k = 0; k < ncp; ++k){
      float u = s[cutpoint[g][k]];
      gout.write((char*)(&u), sizeof(u));
    }
  	gout.close();
  }
  
  int i_sel_id = -1;
  for(int g = 0; g < ngene; ++g){
    R_marg_id[g] = gene_idx[g][marg_id[g]] + 1;
    for(int k = 0; k < (int) sel_id[g].size(); ++k){
      ++i_sel_id;
      R_sel_id[i_sel_id] = gene_idx[g][sel_id[g][k]] + 1;
    }
    int nn = gene_idx[g].size() - sel_id[g].size();
    while(nn){
      ++i_sel_id;
      R_sel_id[i_sel_id] = -1;
      --nn;
    }
  }
  
  int ngap = min(10000, nperm);
  int nblock = nperm / ngap;
  
  #if __NO_PARALLEL__
  srand(seed);
  #endif
  for(int b = 0; b < nblock; ++b){
  	fmat null(ngap, fvec (nsnp, .0f));
  	#if __PARALLEL__
  	drand48_data buf;
  	#endif
  	// compute null score
  	#if __PARALLEL__
  	#pragma omp parallel num_threads(nthread) private(buf)
  	{
  		srand48_r(seed + b * nthread + omp_get_thread_num(), &buf);
  		#pragma omp for
  	#endif
	  	for(int i = 0; i < ngap; ++i){
	  		fvec rn;
	  		#if __PARALLEL__
	  		rnorm(buf, nsnp, rn);
	  		#else
	  		rnorm(nsnp, rn);
	  		#endif
	  		for(int j = 0; j < nsnp; ++j){
	  			null[i][j] = .0f;
	  			for(int k = 0; k < nsnp; ++k){
	  				null[i][j] += rn[k] * U[k][j];
	  			}
	  		}
	  	}
	  #if __PARALLEL__
	  }
	  #endif
	  
	  // write null statistics to local files (per gene)
	  #if __PARALLEL__
	  #pragma omp parallel num_threads(min(nthread, ngene))
	  {
	  	#pragma omp for
	  #endif
	  	for(int g = 0; g < ngene; ++g){
	  		ofstream gout;
	  		gout.open(gene_out[g].c_str(), ios::out | ios::binary | ios::app);
	  		if(!gout){
	  			error("Fail to write null statistics to file");
	  		}
	  		
	  		fmat Sigma;
	  		extract_cov(Sigma, V, gene_idx[g]);
	  		
        int ncp = cutpoint[g].size();
        int mc = cutpoint[g][ncp - 1];
	  		for(int i = 0; i < ngap; ++i){
	  			fvec S;
	  			extract_score(S, null[i], gene_idx[g]);
	  			
	  			fvec s;
          ivec sel_id;
          int marg_id;
          
          if(method == 1){
            search1(s, sel_id, marg_id, S, Sigma, mc);
          }else{
            search2(s, sel_id, marg_id, S, Sigma, mc);
          }
          
          for(int k = 0; k < ncp; ++k){
            float u = s[cutpoint[g][k]];
            gout.write((char*)(&u), sizeof(u));
          }
	  		}
	  		gout.close();
	  	}
	  #if __PARALLEL__
	  }
	  #endif
	  //fmat().swap(null);
  }
  
  // read null statistics (per gene)
  int irk = -1;
  for(int g = 0; g < ngene; ++g){
  	int ncp = cutpoint[g].size();
  	
  	vector<VecStat> stat(ncp, VecStat (nperm + 1, STAT0));
    fstream gin(gene_out[g].c_str(), ios::in | ios::binary);
    for(int i = 0; i < nperm + 1; ++i){
  		for(int j = 0; j < ncp; ++j){
  			float s = .0f;
  			gin.read((char*)(&s), sizeof(s));
  			stat[j][i].stat = s;
  			stat[j][i].id = i;
  		}
  	}
    gin.close();
    
    if(remove(gene_out[g].c_str())){
      error("Cannot delete gene output file");
    }
  	
  	imat arr_rank(ncp, ivec (nperm + 1, 0));
  	#if __PARALLEL__
  	#pragma omp parallel num_threads(min(ncp, nthread))
  	{
      #pragma omp for
    #endif
  		for(int j = 0; j < ncp; ++j){
  			sort(stat[j].begin(), stat[j].end(), descending);
  			for(int i = 0; i < nperm + 1; ++i){
  				int id = stat[j][i].id;
  				arr_rank[j][id] = i;
  			}
  			VecStat().swap(stat[j]);
  		}
  	#if __PARALLEL__
  	}
  	#endif
  	
  	vector<VecStat>().swap(stat);
    
    VecMinp gene_min_p (nperm + 1, MINP0);
    ivec subsum(nthread, 0);
    ivec subtie(nthread, 0);
    int m = nperm + 1;
    for(int j = 0; j < ncp; ++j){
      ++irk;
      R_arr_rank[irk] = arr_rank[j][0];
      if(arr_rank[j][0] < m){
        m = arr_rank[j][0];
      }
    }
    gene_min_p[0].stat = m;
    gene_min_p[0].id = 0;
    
    #if __PARALLEL__
    #pragma omp parallel num_threads(nthread)
    {
      #pragma omp for
    #endif
      for(int i = 1; i < nperm + 1; ++i){
      	#if __PARALLEL__
        int tid = omp_get_thread_num();
        #else
        int tid = 0;
        #endif
        int m = nperm + 1;
        for(int j = 0; j < ncp; ++j){
          if(arr_rank[j][i] < m){
            m = arr_rank[j][i];
          }
        }
        gene_min_p[i].stat = m;
        gene_min_p[i].id = i;
        if(gene_min_p[i].stat < gene_min_p[0].stat){
          subsum[tid] += 1;
        }else if(gene_min_p[i].stat == gene_min_p[0].stat){
          subtie[tid] += 1;
        }else{
          ;
        }
      }
    #if __PARALLEL__
    }
    #endif
    
    R_gene_pval[g] = 1.0;
    int rep = 0;
    for(int t = 0; t < nthread; ++t){
      R_gene_pval[g] += subsum[t];
      rep += subtie[t];
    }
    R_gene_pval[g] += rep / 2.0;
    R_gene_pval[g] /= nperm + 1;
    
    sort(gene_min_p.begin(), gene_min_p.end());
    
    int *gene_p_stat = new int[nperm + 1];
    if(gene_p_stat == NULL){
    	error("Out of memory in artp2_chr");
    }
    
    #if __PARALLEL__
    #pragma omp parallel num_threads(nthread)
    {
    	#pragma omp for
    #endif
	    for(int i = 0; i < nperm + 1; ++i){
	    	int id = gene_min_p[i].id;
	    	gene_p_stat[id] = i;
	    }
	  #if __PARALLEL__
	  }
	  #endif
    VecMinp().swap(gene_min_p);
    
    fstream gout(gene_out[g].c_str(), ios::out | ios::binary);
    if(!gout){
      error("Fail to write gene statistics to file");
    }
    
    gout.write((char*)(gene_p_stat), sizeof(int) * (nperm + 1));
    delete[] gene_p_stat;
    gout.close();
    
  }
  
  delete[] file_prefix;
  
}

void artp2(char **R_file_prefix, int *R_nperm, int *R_nthread, 
int *R_ngene, int *R_group_id, int *R_gene_id, 
int *R_pathway_cutpoint, int *R_ncp, 
double *R_pathway_pval, int *R_arr_rank, double *R_gene_pval){
  
  int len_file_prefix = strlen(*R_file_prefix);
  char *file_prefix = new char[len_file_prefix + 1];
  file_prefix[0] = '\0';
  strcat(file_prefix, *R_file_prefix);
  
  int nperm = *R_nperm;
  int nthread = *R_nthread;
  int ngene = *R_ngene;
  int ncp = *R_ncp;
  
  #if __NO_PARALLEL__
  assert(nthread == 1);
  #endif
  
  ivec group_id;
  load_group_id(R_group_id, group_id, ngene);
  
  ivec gene_id;
  load_gene_id(R_gene_id, gene_id, ngene);
  
  ivec pathway_cutpoint;
  load_pathway_cutpoint(R_pathway_cutpoint, pathway_cutpoint, ncp);
  int max_cutpoint = pathway_cutpoint[ncp - 1];
  
  string fprefix (file_prefix);
  svec gene_out (ngene, fprefix);
  vector<shared_ptr<fstream> > gin;
  for(int g = 0; g < ngene; ++g){
    ostringstream cid;
    cid << group_id[g];
    ostringstream gid;
    gid << gene_id[g];
    gene_out[g] = gene_out[g] + string(".CID.") + cid.str() + string(".GID.") + gid.str() + string(".bin");
    gin.push_back(make_shared<fstream> (gene_out[g].c_str(), ios::in | ios::binary));
  }
  
  
  vector<VecStat> stat(ncp, VecStat (nperm + 1, STAT0));
  int maxlines = 50000;
  int rem = nperm + 1;
  int inp = -1;
  while(rem > 0){
    int nlines = (rem > maxlines) ? maxlines : rem;
  	
  	fmat s(nlines, fvec(ngene, .0f));
  	#if __PARALLEL__
  	#pragma omp parallel num_threads(nthread)
  	{
  		#pragma omp for
  	#endif
	  	for(int g = 0; g < ngene; ++g){
	  		int *buffer = new int[nlines];
	  		if(buffer == NULL){
	  			error("Out of memory in artp2");
	  		}
	  		
	  		(*(gin[g])).read((char*) buffer, nlines * sizeof(int));
	  		if(inp == -1){
	  			R_gene_pval[g] = (buffer[0] + 1.0) / (nperm + 1);
	  		}
	  		
	  		for(int nl = 0; nl < nlines; ++nl){
	  			s[nl][g] = (float) log((buffer[nl] + 1.0) / (nperm + 1));
	  		}
	  		delete[] buffer;
	  	}
	  	
	  	#if __PARALLEL__
  		#pragma omp for
  		#endif
	  	for(int nl = 0; nl < nlines; ++nl){
	  		if(ngene > 1){
	  			sort(s[nl].begin(), s[nl].end());
	  			for(int g = 1; g <= max_cutpoint; ++g){
	  				s[nl][g] += s[nl][g - 1];
	  			}
	  		}
				
				#if __PARALLEL__
				for(int k = 0; k < ncp; ++k){
					float u = -s[nl][pathway_cutpoint[k]];
					stat[k][inp + nl + 1].stat = u;
					stat[k][inp + nl + 1].id = inp + nl + 1;
				}
				#else
				++inp;
				for(int k = 0; k < ncp; ++k){
					float u = -s[nl][pathway_cutpoint[k]];
					stat[k][inp].stat = u;
					stat[k][inp].id = inp;
				}
	  		#endif
				
				fvec().swap(s[nl]);
	  	}
	  #if __PARALLEL__
	  }
	  #endif
	  
	  #if __PARALLEL__
	  inp += nlines;
	  #endif
  	
  	fmat().swap(s);
  	
  	rem -= nlines;
  }
  
  for(int g = 0; g < ngene; ++g){
  	(*(gin[g])).close();
  	if(remove(gene_out[g].c_str())){
  		error("Cannot delete gene output file");
  	}
  }
  
  imat arr_rank(ncp, ivec (nperm + 1, 0));
  #if __PARALLEL__
  #pragma omp parallel num_threads(min(nthread, ncp))
  {
    #pragma omp for
  #endif
    for(int k = 0; k < ncp; ++k){
      sort(stat[k].begin(), stat[k].end(), descending);
      for(int i = 0; i < nperm + 1; ++i){
        int id = stat[k][i].id;
        arr_rank[k][id] = i;
      }
      VecStat().swap(stat[k]);
    }
  #if __PARALLEL__
  }
  #endif
  vector<VecStat>().swap(stat);
  
  ivec pathway_min_p(nperm + 1, -1);
  ivec subsum(nthread, 0);
  ivec subtie(nthread, 0);
  int m = nperm + 1;
  for(int k = 0; k < ncp; ++k){
    R_arr_rank[k] = arr_rank[k][0];
    if(arr_rank[k][0] < m){
      m = arr_rank[k][0];
    }
  }
  pathway_min_p[0] = m;
  
  #if __PARALLEL__
  #pragma omp parallel num_threads(nthread)
  {
    #pragma omp for
  #endif
    for(int i = 1; i < nperm + 1; ++i){
    	#if __PARALLEL__
      int tid = omp_get_thread_num();
      #else
      int tid = 0;
      #endif
      int m = nperm + 1;
      for(int k = 0; k < ncp; ++k){
        if(arr_rank[k][i] < m){
          m = arr_rank[k][i];
        }
      }
      pathway_min_p[i] = m;
      if(pathway_min_p[i] < pathway_min_p[0]){
        subsum[tid] += 1;
      }else if(pathway_min_p[i] == pathway_min_p[0]){
        subtie[tid] += 1;
      }else{
        ;
      }
    }
  #if __PARALLEL__
  }
  #endif
  
  *R_pathway_pval = 1.0;
  int rep = 0;
  for(int t = 0; t < nthread; ++t){
    *R_pathway_pval += subsum[t];
    rep += subtie[t];
  }
  *R_pathway_pval += rep / 2.0;
  *R_pathway_pval /= nperm + 1;
  
  delete[] file_prefix;
  
  
}

}






