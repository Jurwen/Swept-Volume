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
const int timer_amount = 9;

/// The labels for timing stats.
const std::array<std::string, timer_amount> time_label = {"temporal splits",
    "spatial splits",
    "create vertex column",
    "add connectivity",
    "extrude tet column",
    "evaluate tet column",
    "first function",
    "second function",
    "refinement criteria"
};

/// the enum for the timing labels.
enum timeProfileName{
    time_splits,
    space_splits,
    init_vert_col,
    add_connectivity,
    extrude_tet_col,
    eval_tet_col,
    first_func,
    second_func,
    ref_crit
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
            case init_vert_col:
                m_timeProfile[2] += ms;
                break;
            case add_connectivity:
                m_timeProfile[3] += ms;
                break;
            case extrude_tet_col:
                m_timeProfile[4] += ms;
                break;
            case eval_tet_col:
                m_timeProfile[5] += ms;
                break;
            case first_func:
                m_timeProfile[6] += ms;
                break;
            case second_func:
                m_timeProfile[7] += ms;
                break;
            case ref_crit:
                m_timeProfile[8] += ms;
                break;
            default:
                throw std::runtime_error("no matching time profile identifier");
        }
        m_Func(m_timeProfile);
    }
    
    
private:
    timeProfileName m_Name;
    Fn m_Func;
    std::array<double, timer_amount> m_timeProfile = {0,0,0,0,0,0,0,0};
    std::chrono::time_point<std::chrono::high_resolution_clock> starterTime;
};

#endif /* timer_h */
