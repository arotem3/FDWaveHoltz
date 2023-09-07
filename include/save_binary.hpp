#ifndef WH_SAVE_BINARY_HPP
#define WH_SAVE_BINARY_HPP

#include <iostream>
#include <fstream>

// save array of length n to file fname in raw binary format.
void save_binary(const double * x, int n, const std::string& fname);

#endif