/* this code is used to test C++ feature */

#include "util.h"

struct STAT1{
  
  float stat;
  int id;
  
  bool operator < (const STAT1& st) const{
    return (stat < st.stat);
  }
  
  STAT1(float st, int i){
    stat = st;
    id = i;
  }
  
};

STAT1 STAT10(.0f, -1);

bool decreasingf1(float a, float b){
  return a > b;
}

bool decreasing1(const STAT1 &a, const STAT1 &b){
  return a.stat > b.stat;
}

bool increasingf1(float a, float b){
  return a < b;
}

bool increasing1(const STAT1 &a, const STAT1 &b){
  return a.stat < b.stat;
}

extern "C" {

void test(){
  
  float a[] = {1.7, 2.6, 7.9, 11.3, 5.6, 4.3, 1.76};
  vector<float> p (a, a + sizeof(a) / sizeof(float));
  
  cout << "original p-value vector" << endl;
  for(int i = 0; i < p.size(); ++i){
    cout << p[i] << " ";
  }
  cout << endl;
  
  vector<STAT1> stat(p.size(), STAT10);
  cout << "-log(p-value)" << endl;
  for(int i = 0; i < p.size(); ++i){
    stat[i].stat = -pchisq(p[i], 1, false, true);
    stat[i].id = i;
    cout <<  stat[i].stat << " (" << i << ") ";
  }
  cout << endl;
  
  cout << "max p-value = " << *max_element(p.begin(), p.end()) << endl;
  
  vector<float> p0 = p;
  cout << "sorted p-value" << endl;
  sort(p.begin(), p.end());
  for(int i = 0; i < p.size(); ++i){
    cout << p[i] << " ";
  }
  cout << endl;
  
  p = p0;
  cout << "sorted p-value" << endl;
  sort(p.begin(), p.end(), increasingf1);
  for(int i = 0; i < p.size(); ++i){
    cout << p[i] << " ";
  }
  cout << endl;
  
  vector<STAT1> stat0 = stat;
  cout << "sorted stat" << endl;
  sort(stat.begin(), stat.end());
  for(int i = 0; i < stat.size(); ++i){
    cout << stat[i].stat << " (" << stat[i].id << ") ";
  }
  cout << endl;
  
  stat = stat0;
  cout << "sorted stat" << endl;
  sort(stat.begin(), stat.end(), decreasing1);
  for(int i = 0; i < stat.size(); ++i){
    cout << stat[i].stat << " (" << stat[i].id << ") ";
  }
  cout << endl;
  
  stat = stat0;
  cout << "sorted stat" << endl;
  sort(stat.begin(), stat.end());
  for(int i = 0; i < stat.size(); ++i){
    cout << stat[i].stat << " (" << stat[i].id << ") ";
  }
  cout << endl;
  
}

}


