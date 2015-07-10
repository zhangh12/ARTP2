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

void load_sigma2(double *R_sigma2, fvec &sigma2, int nsnp){
  
  sigma2 = fvec (nsnp, .0f);
	for(int i = 0; i < nsnp; ++i){
		sigma2[i] = R_sigma2[i];
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


#if __PARALLEL__
extern "C" {

void artp3_chr(char **R_file_prefix, int *R_nperm, int *R_seed, 
int *R_nthread, int *R_nsnp, int *R_ngene, 
double *R_vU, double *R_score0, double *R_sigma2, 
int *R_vgene_idx, int *R_gene_start, int *R_gene_end, 
int *R_vgene_cutpoint, 
int *R_gene_cutpoint_start, int *R_gene_cutpoint_end, 
double *R_gene_pval, int *R_arr_rank){
  
  int len_file_prefix = strlen(*R_file_prefix);
  char *file_prefix = new char[len_file_prefix + 1];
  file_prefix[0] = '\0';
  strcat(file_prefix, *R_file_prefix);
  
  int nperm = *R_nperm;
  int seed = *R_seed;
  int nthread = *R_nthread;
  int nsnp = *R_nsnp;
  int ngene = *R_ngene;
  
  fvec score0;
  fvec sigma2;
  fmat U;
  load_score0(R_score0, score0, nsnp);
  load_sigma2(R_sigma2, sigma2, nsnp);
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
  for(int g = 0; g < ngene; ++g){
  	fstream gout(gene_out[g].c_str(), ios::out | ios::binary);
  	if(!gout){
  		error("Fail to write null statistics to file");
  	}
  	int ns = gene_idx[g].size();
    int ncp = cutpoint[g].size();
    int max_cutpoint = cutpoint[g][ncp - 1];
    fvec s (ns, .0f);
  	for(int j = 0; j < ns; ++j){
  		s[j] = score0[gene_idx[g][j]];
      s[j] = pchisq(s[j] * s[j] / sigma2[gene_idx[g][j]], 1, false, true);
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
  
  int ngap = min(10000, nperm);
  int nblock = nperm / ngap;
  
  for(int b = 0; b < nblock; ++b){
  	fmat null(ngap, fvec (nsnp, .0f));
  	drand48_data buf;
  	// compute null statistics
  	#pragma omp parallel num_threads(nthread) private(buf)
  	{
  		srand48_r(seed + b * nthread + omp_get_thread_num(), &buf);
  		#pragma omp for
	  	for(int i = 0; i < ngap; ++i){
	  		fvec rn;
	  		rnorm(buf, nsnp, rn);
	  		for(int j = 0; j < nsnp; ++j){
	  			null[i][j] = .0f;
	  			for(int k = 0; k < nsnp; ++k){
	  				null[i][j] += rn[k] * U[k][j];
	  			}
	  			null[i][j] = null[i][j] * null[i][j] / sigma2[j];
	  			null[i][j] = pchisq(null[i][j], 1, false, true);
	  		}
	  	}
	  }
	  
	  // write null statistics to local files (per gene)
	  #pragma omp parallel num_threads(min(nthread, ngene))
	  {
	  	#pragma omp for
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
	  }
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
  	
  	imat arr_rank(ncp, ivec (nperm + 1, 0));
  	#pragma omp parallel num_threads(min(ncp, nthread))
  	{
      #pragma omp for
  		for(int j = 0; j < ncp; ++j){
  			sort(stat[j].begin(), stat[j].end(), descending);
  			for(int i = 0; i < nperm + 1; ++i){
  				int id = stat[j][i].id;
  				arr_rank[j][id] = i;
  			}
  		}
  	}
  	
  	vector<VecStat>().swap(stat);
    
    ivec gene_min_p (nperm + 1, -1);
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
    gene_min_p[0] = m;
    
    #pragma omp parallel num_threads(nthread)
    {
      #pragma omp for
      for(int i = 1; i < nperm + 1; ++i){
        int tid = omp_get_thread_num();
        int m = nperm + 1;
        for(int j = 0; j < ncp; ++j){
          if(arr_rank[j][i] < m){
            m = arr_rank[j][i];
          }
        }
        gene_min_p[i] = m;
        if(gene_min_p[i] < gene_min_p[0]){
          subsum[tid] += 1;
        }else if(gene_min_p[i] == gene_min_p[0]){
          subtie[tid] += 1;
        }else{
          ;
        }
      }
    }
    
    R_gene_pval[g] = 1.0;
    int rep = 0;
    for(int t = 0; t < nthread; ++t){
      R_gene_pval[g] += subsum[t];
      rep += subtie[t];
    }
    R_gene_pval[g] += rep / 2.0;
    R_gene_pval[g] /= nperm + 1;
    
    ofstream gout;
    gout.open(gene_out[g].c_str(), ios::out | ios::binary);
    if(!gout){
      error("Fail to write gene statistics to file");
    }
    for(int i = 0; i < nperm + 1; ++i){
      gout.write((char*)(&(gene_min_p[i])), sizeof(gene_min_p[i]));
    }
    gout.close();
    
  }
  
  
  delete[] file_prefix;
  
}

void artp3(char **R_file_prefix, int *R_nperm, int *R_nthread, 
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
  
  ivec group_id;
  load_group_id(R_group_id, group_id, ngene);
  
  ivec gene_id;
  load_gene_id(R_gene_id, gene_id, ngene);
  
  ivec pathway_cutpoint;
  load_pathway_cutpoint(R_pathway_cutpoint, pathway_cutpoint, ncp);
  int max_cutpoint = pathway_cutpoint[ncp - 1];
  
  string fprefix (file_prefix);
  svec gene_out (ngene, fprefix);
  for(int g = 0; g < ngene; ++g){
    ostringstream cid;
    cid << group_id[g];
    ostringstream gid;
    gid << gene_id[g];
    gene_out[g] = gene_out[g] + string(".CID.") + cid.str() + string(".GID.") + gid.str() + string(".bin");
  }
  
  vector<VecMinp> gene_stat (ngene, VecMinp(nperm + 1, MINP0));
  imat gene_p_stat (ngene, ivec(nperm + 1, -1));
  #pragma omp parallel num_threads(min(nthread, ngene))
  {
    #pragma omp for
    for(int g = 0; g < ngene; ++g){
      fstream gin(gene_out[g].c_str(), ios::in | ios::binary);
      //ivec zero_loc;
      for(int i = 0; i < nperm + 1; ++i){
        gin.read((char*)(&(gene_stat[g][i].stat)), sizeof(gene_stat[g][i].stat));
        gene_stat[g][i].id = i;
      }
      gin.close();
      if(remove(gene_out[g].c_str())){
        error("Cannot delete gene output file");
      }
      sort(gene_stat[g].begin(), gene_stat[g].end());
      for(int i = 0; i < nperm + 1; ++i){
        int id = gene_stat[g][i].id;
        gene_p_stat[g][id] = i;
      }
      VecMinp().swap(gene_stat[g]);
    }
  }
  
  vector<VecMinp>().swap(gene_stat);
  
  for(int g = 0; g < ngene; ++g){
    R_gene_pval[g] = (gene_p_stat[g][0] + 1.0) / (nperm + 1);
  }
  
  vector<VecStat> stat(ncp, VecStat (nperm + 1, STAT0));
  #pragma omp parallel num_threads(nthread)
  {
    #pragma omp for
    for(int i = 0; i < nperm + 1; ++i){
      fvec s(ngene, .0f);
      for(int g = 0; g < ngene; ++g){
        s[g] = (float) log((gene_p_stat[g][i] + 1.0) / (nperm + 1));
      }
      if(ngene > 1){
        sort(s.begin(), s.end());
        for(int g = 1; g <= max_cutpoint; ++g){
          s[g] += s[g - 1];
        }
      }
      
      for(int k = 0; k < ncp; ++k){
        float u = -s[pathway_cutpoint[k]];
        stat[k][i].stat = u;
        stat[k][i].id = i;
      }
    }
  }
  
  imat().swap(gene_p_stat);
  
  imat arr_rank(ncp, ivec (nperm + 1, 0));
  #pragma omp parallel num_threads(min(nthread, ncp))
  {
    #pragma omp for
    for(int k = 0; k < ncp; ++k){
      sort(stat[k].begin(), stat[k].end(), descending);
      for(int i = 0; i < nperm + 1; ++i){
        int id = stat[k][i].id;
        arr_rank[k][id] = i;
      }
    }
  }
  
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
  
  #pragma omp parallel num_threads(nthread)
  {
    #pragma omp for
    for(int i = 1; i < nperm + 1; ++i){
      int tid = omp_get_thread_num();
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
  }
  
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
#else
extern "C" {

void artp3_chr(char **R_file_prefix, int *R_nperm, int *R_seed, 
int *R_nthread, int *R_nsnp, int *R_ngene, 
double *R_vU, double *R_score0, double *R_sigma2, 
int *R_vgene_idx, int *R_gene_start, int *R_gene_end, 
int *R_vgene_cutpoint, 
int *R_gene_cutpoint_start, int *R_gene_cutpoint_end, 
double *R_gene_pval, int *R_arr_rank){
  
  int len_file_prefix = strlen(*R_file_prefix);
  char *file_prefix = new char[len_file_prefix + 1];
  file_prefix[0] = '\0';
  strcat(file_prefix, *R_file_prefix);
  
  int nperm = *R_nperm;
  int seed = *R_seed;
  int nthread = *R_nthread;
  int nsnp = *R_nsnp;
  int ngene = *R_ngene;
  
  assert(nthread == 1);
  
  fvec score0;
  fvec sigma2;
  fmat U;
  load_score0(R_score0, score0, nsnp);
  load_sigma2(R_sigma2, sigma2, nsnp);
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
  for(int g = 0; g < ngene; ++g){
    fstream gout(gene_out[g].c_str(), ios::out | ios::binary);
  	if(!gout){
  		error("Fail to write null statistics to file");
  	}
  	int ns = gene_idx[g].size();
    int ncp = cutpoint[g].size();
    int max_cutpoint = cutpoint[g][ncp - 1];
    fvec s (ns, .0f);
  	for(int j = 0; j < ns; ++j){
  		s[j] = score0[gene_idx[g][j]];
      s[j] = pchisq(s[j] * s[j] / sigma2[gene_idx[g][j]], 1, false, true);
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
  
  int ngap = min(10000, nperm);
  int nblock = nperm / ngap;
  
  srand(seed);
  for(int b = 0; b < nblock; ++b){
  	fmat null(ngap, fvec (nsnp, .0f));
  	// compute null statistics
  	for(int i = 0; i < ngap; ++i){
  		fvec rn;
  		rnorm(nsnp, rn);
  		for(int j = 0; j < nsnp; ++j){
  			null[i][j] = .0f;
  			for(int k = 0; k < nsnp; ++k){
  				null[i][j] += rn[k] * U[k][j];
  			}
  			null[i][j] = null[i][j] * null[i][j] / sigma2[j];
  			null[i][j] = pchisq(null[i][j], 1, false, true);
  		}
  	}
	  
	  // write null statistics to local files (per gene)
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
  	
  	imat arr_rank(ncp, ivec (nperm + 1, 0));
		for(int j = 0; j < ncp; ++j){
			sort(stat[j].begin(), stat[j].end(), descending);
			for(int i = 0; i < nperm + 1; ++i){
				int id = stat[j][i].id;
				arr_rank[j][id] = i;
			}
		}
  	
  	vector<VecStat>().swap(stat);
    
    ivec gene_min_p (nperm + 1, -1);
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
    gene_min_p[0] = m;
    
    for(int i = 1; i < nperm + 1; ++i){
      int tid = 0;
      int m = nperm + 1;
      for(int j = 0; j < ncp; ++j){
        if(arr_rank[j][i] < m){
          m = arr_rank[j][i];
        }
      }
      gene_min_p[i] = m;
      if(gene_min_p[i] < gene_min_p[0]){
        subsum[tid] += 1;
      }else if(gene_min_p[i] == gene_min_p[0]){
        subtie[tid] += 1;
      }else{
        ;
      }
    }
    
    R_gene_pval[g] = 1.0;
    int rep = 0;
    for(int t = 0; t < nthread; ++t){
      R_gene_pval[g] += subsum[t];
      rep += subtie[t];
    }
    R_gene_pval[g] += rep / 2.0;
    R_gene_pval[g] /= nperm + 1;
    
    ofstream gout;
    gout.open(gene_out[g].c_str(), ios::out | ios::binary);
    if(!gout){
      error("Fail to write gene statistics to file");
    }
    for(int i = 0; i < nperm + 1; ++i){
      gout.write((char*)(&(gene_min_p[i])), sizeof(gene_min_p[i]));
    }
    gout.close();
    
  }
  
  
  delete[] file_prefix;
  
}

void artp3(char **R_file_prefix, int *R_nperm, int *R_nthread, 
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
  
  assert(nthread == 1);
  
  ivec group_id;
  load_group_id(R_group_id, group_id, ngene);
  
  ivec gene_id;
  load_gene_id(R_gene_id, gene_id, ngene);
  
  ivec pathway_cutpoint;
  load_pathway_cutpoint(R_pathway_cutpoint, pathway_cutpoint, ncp);
  int max_cutpoint = pathway_cutpoint[ncp - 1];
  
  string fprefix (file_prefix);
  svec gene_out (ngene, fprefix);
  for(int g = 0; g < ngene; ++g){
    ostringstream cid;
    cid << group_id[g];
    ostringstream gid;
    gid << gene_id[g];
    gene_out[g] = gene_out[g] + string(".CID.") + cid.str() + string(".GID.") + gid.str() + string(".bin");
  }
  
  vector<VecMinp> gene_stat (ngene, VecMinp(nperm + 1, MINP0));
  imat gene_p_stat (ngene, ivec(nperm + 1, -1));
  
  for(int g = 0; g < ngene; ++g){
    fstream gin(gene_out[g].c_str(), ios::in | ios::binary);
    //ivec zero_loc;
    for(int i = 0; i < nperm + 1; ++i){
      gin.read((char*)(&(gene_stat[g][i].stat)), sizeof(gene_stat[g][i].stat));
      gene_stat[g][i].id = i;
    }
    gin.close();
    if(remove(gene_out[g].c_str())){
      error("Cannot delete gene output file");
    }
    sort(gene_stat[g].begin(), gene_stat[g].end());
    for(int i = 0; i < nperm + 1; ++i){
      int id = gene_stat[g][i].id;
      gene_p_stat[g][id] = i;
    }
    VecMinp().swap(gene_stat[g]);
  }
  
  vector<VecMinp>().swap(gene_stat);
  
  for(int g = 0; g < ngene; ++g){
    R_gene_pval[g] = (gene_p_stat[g][0] + 1.0) / (nperm + 1);
  }
  
  vector<VecStat> stat(ncp, VecStat (nperm + 1, STAT0));
  for(int i = 0; i < nperm + 1; ++i){
    fvec s(ngene, .0f);
    for(int g = 0; g < ngene; ++g){
      s[g] = (float) log((gene_p_stat[g][i] + 1.0) / (nperm + 1));
    }
    if(ngene > 1){
      sort(s.begin(), s.end());
      for(int g = 1; g <= max_cutpoint; ++g){
        s[g] += s[g - 1];
      }
    }
    
    for(int k = 0; k < ncp; ++k){
      float u = -s[pathway_cutpoint[k]];
      stat[k][i].stat = u;
      stat[k][i].id = i;
    }
  }
  
  imat().swap(gene_p_stat);
  
  imat arr_rank(ncp, ivec (nperm + 1, 0));
  for(int k = 0; k < ncp; ++k){
    sort(stat[k].begin(), stat[k].end(), descending);
    for(int i = 0; i < nperm + 1; ++i){
      int id = stat[k][i].id;
      arr_rank[k][id] = i;
    }
  }
  
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
  
  for(int i = 1; i < nperm + 1; ++i){
    int tid = 0;
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
#endif

