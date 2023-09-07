#ifndef PROGRESS_BAR_HPP
#define PROGRESS_BAR_HPP

#include <string>

class ProgressBar
{
    private:
    std::string bar;
    int k;
    int step;
    int j;
    
    public:
    ProgressBar(int total, int len=30) : bar(len, ' '), k{0}, step{total/len}, j{0} {}

    const std::string& get_progress() const
    {
        return bar;
    }

    ProgressBar& operator++()
    {
        j = (j+1)%step;
        if (j == 0 && (size_t)k < bar.size())
        {
            bar.at(k) = '#';
            k++;
        }
        return *this;
    }
};


#endif