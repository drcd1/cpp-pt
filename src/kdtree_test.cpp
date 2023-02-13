#include <algorithms/kdtree.h>
#include <cstdlib>
#include <vector>

using namespace cpppt;

struct kdtreept{
    Vec3 pt;
};

struct mark{
    int visited = 0;
    int inrange = 0;
    int bf_inrange = 0;
    float distance;

};

#define RANDF() float(rand())/float(RAND_MAX)

int main(){
    KDTree<kdtreept> t;
    std::vector<mark> marks(100);

    for(int i = 0; i<100; i++){
        kdtreept p;
        p.pt = Vec3(RANDF(),RANDF(),RANDF());
        t.add(p);
    }

    t.build();

    Vec3 pt = Vec3(RANDF(),RANDF(),RANDF());
    float radius = 0.2;
    t.run_kernel(pt,radius,[&](int other){
        marks.at(other).visited = 1;
        if(lensqr(pt-t.points.at(other).pt)<radius*radius){
            marks.at(other).inrange = 1;
        }
    });

    for(int i = 0; i<100; i++){
        Vec3 opt = t.points.at(i).pt;
        if(lensqr(pt-opt)<radius*radius){
            marks.at(i).bf_inrange = 1;
        }
        marks.at(i).distance = length(pt-opt);

        std::cout<<(marks.at(i).visited?std::string("yes | "):std::string("no  | "))<<
            (marks.at(i).inrange?std::string("yes | "):std::string("no  | "))<<
            (marks.at(i).bf_inrange?std::string("yes | "):std::string("no  | "))<<
            (marks.at(i).distance)<<std::endl;
    }



}