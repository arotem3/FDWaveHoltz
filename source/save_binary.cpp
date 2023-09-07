#include "save_binary.hpp"

void save_binary(const double * x, int n, const std::string& fname)
{
    std::ofstream out(fname, std::ios::out | std::ios::binary);
    out.write(reinterpret_cast<const char*>(x), n * sizeof(double));
    out.close();
}