#ifndef LOADMTX_H__
#define LOADMTX_H__

#include <vector>
#include <string>

#include "mtxvec.h"

Matrix loadMTXFile(const std::string&);
void writeMTXFile(const std::string&, std::vector<std::vector<double>>);

#endif
