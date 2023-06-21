#ifndef renderer_interface
#define renderer_interface

#include <string>

enum SetingsAceessRule{
    NONE
};

struct SettingProperties{
    SettingsAccessRule rule;
    std::string name;

};


class RendererInterface{

public:
    
    //setting are setby the gui and read in the renderer
    int addSetting(void* setting, Type type, SettngProperties setingProperties){

    }
    //statistics are set in the renderer and read by the gui
    int addStatistic(void* statistic, Type type, StatisticsAccessRule rule){

    }

    //use case: a metropolis renderer records a chain
    //this chain can be played back in the gui
    //this chain can be seen in realtime as the renderer goes



    void GUISettings(){

    }

    void GUIStatistics(){

    }

}







#endif