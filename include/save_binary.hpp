#ifndef WH_SAVE_BINARY_HPP
#define WH_SAVE_BINARY_HPP

#include <iostream>
#include <fstream>

// save array of length n to file fname in raw binary format.
template <typename type>
bool save_binary(const type * x, int n, const std::string& fname)
{
    std::ofstream out(fname, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "error: could not open file " << fname << " for writing\n";
        return false;
    }

    out.write(reinterpret_cast<const char*>(x), n * sizeof(type));
    out.close();

    return true;
}

#endif