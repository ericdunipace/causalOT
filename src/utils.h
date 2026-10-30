#ifndef UTILS_H
#define UTILS_H
#include "causalOT_types.h"

double logSumExp(const vector & x_);
double logSumExp(const matrix & x_);
vector rowLogSumExp(const matrix & x_);
vector colLogSumExp(const matrix & x_);

#endif
