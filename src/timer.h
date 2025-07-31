//
//  timer.h
//  swept_volume
//
//  Created by Yiwen Ju on 12/23/24.
//

#ifndef timer_h
#define timer_h

#include <chrono>
#include <array>
#include <string>

#define time_profile 0

///The current amount in time profiling.
const int timer_amount = 11;

/// The labels for timing stats.
const std::array<std::string, timer_amount> time_label = {"first function",
    "first part set up",
    "first function bezier",
    "non simple poly",
    "compute caps",
    "extract first isosurface",
    "find intersect",
    "first part",
    "second function",
    "refinement criteria",
    "second part"
};

/// the enum for the timing labels.
enum timeProfileName{
    first_func,
    first_part_setup,
    first_bezier,
    non_simple_poly,
    compute_caps,
    extract_first_iso,
    find_intersect,
    first_part,
    second_func,
    ref_crit,
    second_part
};

/// add the timer recorded to the profiling timer.
/// @param[in] profile          The most current time profile.
/// @param[in] timer            The time recorded from this temporary timer.
///
/// @return         The updated time profile.
inline void combine_timer(std::array<double, timer_amount>& profile, std::array<size_t, timer_amount>& profileCount, const size_t& index, const double& ms){
    profile[index] += ms;
    profileCount[index] ++;
}

template<typename Fn>
class Timer
{
public:
    Timer(timeProfileName name,
          Fn&& func
          )
    : m_Name(name), m_Func(func)
    {
        starterTime = std::chrono::high_resolution_clock::now();
    }
    
    ~Timer(){}
    
    void Stop()
    {
        auto stopperTime = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
        auto duration = end - start;
        ms = duration * 1e-3;
        m_Func(m_Name, ms);
    }
    
    
private:
    timeProfileName m_Name;
    Fn m_Func;
    double ms = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> starterTime;
};

#endif /* timer_h */
