#ifndef CPPPT_SAMPLER
#define CPPPT_SAMPLER
#include <random>

namespace cpppt{

class Sampler{
public:
    virtual float sample() = 0;
};

class RandomSampler : public Sampler{
    std::mt19937 generator;
    std::uniform_real_distribution<float> distribution;

public:
    RandomSampler(int seed = 0):generator(seed),distribution(0.0,1.0){
    }
    float sample(){
        return distribution(generator);
    }
};

}

#endif