#include <iostream>
#include <vector>
#include <chrono>
#include <map>
#include <fstream>
#include "type.h"

#define TIME_BEGIN(a) auto time_begin_##a = std::chrono::high_resolution_clock::now()

#define TIME_END(a)   auto time_end_##a = std::chrono::high_resolution_clock::now();\
					  auto elapse_##a = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_##a - time_begin_##a);\
                      RTLOG("[%s time measured : %.5f seconds.]\n", #a, elapse_##a.count() * 1e-9)

namespace rtfmm
{
class Timer
{
public:
    void begin(std::string name)
    {
        name2index[name] = begins.size();
        names.push_back(name);
        begins.push_back(std::chrono::high_resolution_clock::now());
    }
    void end(std::string name)
    {
        auto end = std::chrono::high_resolution_clock::now();
        int index = name2index[name];
        auto elapse = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begins[index]);
        times.push_back(elapse.count() * 1e-9);
    }
    void save(std::string path)
    {
        std::ofstream outputfile(path, std::ios::app);
        for(int i = 0; i < times.size(); i++)
        {
            std::string ss = format("[%s time measured : %f seconds.]", names[i].c_str(), times[i]);
            outputfile << ss << std::endl;
        }
        outputfile << std::endl;
        outputfile.close();
    }
private:
    std::vector<std::chrono::high_resolution_clock::time_point> begins;
    std::vector<float> times;
    std::vector<std::string> names;
    std::map<std::string, int> name2index;
};
}