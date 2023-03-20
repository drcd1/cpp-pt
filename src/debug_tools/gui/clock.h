#ifndef CPPPT_DEBUG_TOOLS_CLOCK
#define CPPPT_DEBUG_TOOLS_CLOCK

#include <chrono>

namespace cpppt{
namespace DebugTools{
class Clock{
    using time_point = std::chrono::time_point<std::chrono::high_resolution_clock>;


    time_point start_t;
    time_point old_t;
    float time;
    float dt;
    bool running;
public:
    Clock(){
        running= false;
    }

    bool isRunning(){
        return running;
    }

    void start(){
        running=true;
        start_t = std::chrono::high_resolution_clock::now();
        old_t = start_t;

    }

    float getDt(){
        return dt;
    }


    void tick(){
        time_point current = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration<float, std::chrono::seconds::period>(current - start_t).count();
        dt = std::chrono::duration<float, std::chrono::seconds::period>(current - old_t).count();
        old_t = current;
    }


};

}


}


#endif