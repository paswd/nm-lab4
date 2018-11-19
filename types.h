#ifndef TYPES_H
#define TYPES_H

#include <limits>
#include <vector>

typedef long double TNum;
typedef std::vector<TNum> TVector;

const TNum EPS = (TNum) std::numeric_limits<TNum>::epsilon() * 10000000.;
const TNum PI = 3.14159265358979323846;

#endif //TYPES_H