#ifndef LOADMTX_H__
#define LOADMTX_H__

#include <vector>
#include <string>

std::vector<std::vector<double>> loadMTXFile(std::string);
void writeMTXFile(std::string, std::vector<std::vector<double>>);

#endif
