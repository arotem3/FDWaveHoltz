#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

class Timer
{
    std::chrono::_V2::system_clock::time_point tic;

    public:
    Timer() : tic(std::chrono::high_resolution_clock::now()) {}

    // seconds since construction
    double elapsed() const
    {
        auto toc = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
        return dur.count() * 1e-9;
    }

    void restart()
    {
        tic = std::chrono::high_resolution_clock::now();
    }
};

#endif