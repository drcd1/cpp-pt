#ifndef CPPPT_SAMPLER
#define CPPPT_SAMPLER
#include <random>
#include <vector>

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

//TODO: do we need the old samples?

class PMCSampler: public Sampler {
    std::vector<float> samples;
    std::vector<float> old_samples;
    RandomSampler* s=nullptr;
    int i = 0;
    float bigMutateProbability = 0.1;
public:

    PMCSampler(){}
    PMCSampler(RandomSampler* s): s(s){}
    
    void set_sampler(RandomSampler* s_a){
        s = s_a;
    }

    float sample(){
        while(i>=samples.size()){
            samples.push_back(s->sample());
        }
        i++;
        return samples.at(i-1);
    }

    void accept(){
        //do nothing
        i=0;
    }

    //todo: move? instead of copy?
    void reject(){
        i=0;
        samples = old_samples;
    }

    void mutate(){
        i=0;
        old_samples = samples;
        float ecs = s->sample();
        if(ecs>bigMutateProbability){
                //small mutate
            for(int i = 0; i<samples.size(); i++){
                float idc; //i don't care
                //todo: modff?
                samples.at(i) = std::modf(samples.at(i) + 0.06*(s->sample()-0.5),&idc);
                samples.at(i) = samples.at(i)<0.0?samples.at(i)+1.0:samples.at(i);
            }
        } else {
            for(int i = 0; i<samples.size(); i++){
                samples.at(i) = s->sample();
            }
        }

    }
    const std::vector<float>& sample_vector(){
        return samples;
    }
    const std::vector<float>& old_sample_vector(){
        return old_samples;
    }

    void set_sample(const std::vector<float>& sample){
        samples = sample;
    }

};
    

class MLTSampler: public Sampler{
    std::vector<float> samples;
    std::vector<float> old_samples;
    RandomSampler s;
    int i = 0;
    float bigMutateProbability = 0.1;
public:
    MLTSampler(int seed): s(seed){}
    float sample(){
        while(i>=samples.size()){
            samples.push_back(s.sample());
        }
        i++;
        return samples.at(i-1);
    }

    void accept(){
        //do nothing
        i=0;
    }

    //todo: move? instead of copy?
    void reject(){
        i=0;
        samples = old_samples;
    }

    void mutate(){
        i=0;
        old_samples = samples;
        float ecs = s.sample();
        if(ecs>bigMutateProbability){
                //small mutate
            for(int i = 0; i<samples.size(); i++){
                float idc; //i don't care
                //todo: modff?
                samples.at(i) = std::modf(samples.at(i) + 0.06*(s.sample()-0.5),&idc);
                samples.at(i) = samples.at(i)<0.0?samples.at(i)+1.0:samples.at(i);
            }
        } else {
            for(int i = 0; i<samples.size(); i++){
                samples.at(i) = s.sample();
            }
        }

    }
    const std::vector<float>& sample_vector(){
        return samples;
    }
    const std::vector<float>& old_sample_vector(){
        return old_samples;
    }

    void set_sample(const std::vector<float>& sample){
        samples = sample;
    }

};

}

#endif