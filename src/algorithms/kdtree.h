#ifndef CPPPT_KDTREE
#define CPPPT_KDTREE
#include <math/math.h>

#include <vector>
#include <algorithm>
namespace cpppt{

template
<typename KDTreePoint>
class KDTree{
private:

    struct Node{
        Node(int left, int right, int id): left(left),right(right),id(id){}

        int left;
        int right;
        int id;
    };
    std::vector<Node> nodes;

    int root;

    float get_axis(Vec3 pt,int axis){
        switch(axis){
            case 0:
                return pt.x;
            case 1:
                return pt.y;
            case 2:
                return pt.z;
            default:
                throw std::runtime_error("Should not happen");
                break;
        }
    }

public:
    KDTree(){

    }

    std::vector<KDTreePoint> points;

    void add(KDTreePoint p){
        points.push_back(p);
        nodes.push_back(Node(-1,-1,points.size()-1));
    }

    void build(){
        std::vector<int> pts(points.size());
        for(int i = 0; i<points.size(); i++){
            pts.at(i) = i;
        }

        root = build_node(pts,0);
    }

    #define gt(a,b) get_axis(points.at(a).pt,axis)>get_axis(points.at(b).pt,axis)
    #define lt(a,b) get_axis(points.at(a).pt,axis)<get_axis(points.at(b).pt,axis)
    #define swap(a,b) int tmp=a; a=b; b=tmp;
    int order_node3(std::vector<int>& pts, int axis){
        //return median
        if(pts.size()==0){
            return -1;
        }
        if(pts.size()==1){
            return pts.at(0);
        }
        if(pts.size()==2){
            if(gt(pts.at(1),pts.at(0))){
                nodes.at(pts.at(1)).left = pts.at(0);
                return pts.at(1);
            } else {
                nodes.at(pts.at(0)).left = pts.at(1);
                return pts.at(0);
            }
        }
        if(gt(pts.at(0),pts.at(1))){
            swap(pts.at(0),pts.at(1));
        }

        if(gt(pts.at(1),pts.at(2))){
            swap(pts.at(1),pts.at(1));
        }
        if(gt(pts.at(0),pts.at(1))){
            swap(pts.at(0),pts.at(1));
        }
        nodes.at(pts.at(1)).left = pts.at(0);
        nodes.at(pts.at(1)).right = pts.at(2);
        return pts.at(1);
    }


    int guess_median(std::vector<int>& pts,int axis){
        if(pts.size()<20){
            std::sort(pts.begin(),pts.end(),[&](int a, int b){
                return (lt(a,b));
            });
            return pts.at(pts.size()/2);
        } else {
            std::vector<int> pt_samples(20);
            int step = pts.size()/20;
            for(int i = 0; i<20; i++){
                //todo: choose randomly?
                pt_samples.at(i) = pts.at(i*step);
            }

            std::sort(pt_samples.begin(),pt_samples.end(),[&](int a, int b){
                return (lt(a,b));
            });
            return pt_samples.at(pt_samples.size()/2);
        }
    }



    int build_node(std::vector<int>& pts, int depth){
        int axis = depth%3;

        if(pts.size()<=3){

            return order_node3(pts,axis);
        }
        else{
            int median = guess_median(pts,axis);
            std::vector<int> left;
            std::vector<int> right;
            for(int i = 0; i<pts.size(); i++){
                if(pts.at(i)==median){
                    continue;
                }
                if(lt(pts.at(i),median)) {
                    left.push_back(pts.at(i));
                } else { //todo: make this random
                    right.push_back(pts.at(i));
                }
            }

            nodes.at(median).left = build_node(left,depth+1);
            nodes.at(median).right = build_node(right,depth+1);
            return median;



        }

    }

    #undef gt
    #undef lt
    #undef swap

    template <typename F>
    void run_kernel(Vec3 pt, float radius, F kernel){
        float r_sqr = radius*radius;
        std::vector<int> node_queue;
        std::vector<int> depths;
        node_queue.push_back(root);
        depths.push_back(0);
        while(!node_queue.empty()){
            int a = node_queue.back();
            node_queue.pop_back();
            int d = depths.back();
            depths.pop_back();

            int axis = d%3;
            //if(dot(points.at(a).pt - pt)<r_sqr){
                kernel(a);
            //}
            if(get_axis(points.at(a).pt,axis)>get_axis(pt,axis)-radius){
                if(nodes.at(a).left>-1){
                    node_queue.push_back(nodes.at(a).left);
                    depths.push_back(d+1);
                }
            }
            if(get_axis(points.at(a).pt,axis)<get_axis(pt,axis)+radius){
                if(nodes.at(a).right>-1){
                    node_queue.push_back(nodes.at(a).right);
                    depths.push_back(d+1);
                }
            }
        }
    }

};

}
#endif