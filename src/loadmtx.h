#ifndef LOADMTX_H__
#define LOADMTX_H__

#include <vector>
#include <string>

#include "mtxvec.h"

Matrix loadMTXFile(const std::string&);
void writeVecToMTXFile(const std::string&, const Vector&);
void writeMatToMTXFile(const std::string&, const Matrix&);

#endif
