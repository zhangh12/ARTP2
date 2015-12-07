#ifndef UTIL_H
#define UTIL_H

#include <Rconfig.h>

#if defined(SUPPORT_OPENMP) && __GNUG__ && defined(__linux__)
#define __PARALLEL__ true
#else
#define __PARALLEL__ false
#endif

#if __PARALLEL__
#define __NO_PARALLEL__ false
#else
#define __NO_PARALLEL__ true
#endif

#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>

#if __PARALLEL__
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

using namespace std;


#endif



