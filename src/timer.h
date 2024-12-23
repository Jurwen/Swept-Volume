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
///The current amount in time profiling.
const int timer_amount = 2;

/// The labels for timing stats.
/// an array of 10 timings: {total time getting the multiple indices, total time,time spent on single function, time spent on double functions, time spent on triple functions time spent on double functions' zero crossing test, time spent on three functions' zero crossing test, total subdivision time, total evaluation time,total splitting time}
const std::array<std::string, timer_amount> time_label = {"temporal splits: ",
    "spatial splits: ",
};

/// the enum for the timing labels.
enum timeProfileName{
    time_splits,
    space_splits
};

/// add the timer recorded to the profiling timer.
/// @param[in] profile          The most current time profile.
/// @param[in] timer            The time recorded from this temporary timer.
///
/// @return         The updated time profile.
std::array<double, timer_amount> combine_timer (const std::array<double, timer_amount> &profile, const std::array<double, timer_amount> &timer);

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
        double ms = duration * 0.001;
        switch (m_Name){
            case time_splits:
                m_timeProfile[0] += ms;
                break;
            case space_splits:
                m_timeProfile[1] += ms;
                break;
            default:
                throw std::runtime_error("no matching time profile identifier");
        }
        m_Func(m_timeProfile);
    }
    
    
private:
    timeProfileName m_Name;
    Fn m_Func;
    std::array<double, timer_amount> m_timeProfile = {0,0};
    std::chrono::time_point<std::chrono::high_resolution_clock> starterTime;
};

#endif /* timer_h */
