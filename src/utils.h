#ifndef UTILS_H
#define UTILS_H
#include "causalOT_types.h"

double logSumExp(vector & x_);
double logSumExp(matrix & x_);
vector rowLogSumExp(matrix & x_);
vector colLogSumExp(matrix & x_);

#endif