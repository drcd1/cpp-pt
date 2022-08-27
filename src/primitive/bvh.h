#ifndef CPPPT_BVH
#define CPPPT_BVH

#include <primitive/ray.h>
#include <primitive/primitive.h>
#include <primitive/simple_group.h>
#include <primitive/aabb.h>
#include <memory>


namespace cpppt{

struct BVHNode: public Primitive{
    std::unique_ptr<Primitive> right;
    std::unique_ptr<Primitive> left;
    AABB aabb;

    BVHNode(std::unique_ptr<Primitive>& right, std::unique_ptr<Primitive>& left, AABB& aabb): right(std::move(right)),left(std::move(left)),aabb(aabb){

    }

    virtual bool intersect(Ray& r, Intersection* is) const {
        if(!aabb.intersect_any(r))
            return false;

        #ifdef RAY_STATISTICS
        r.tests+=1;
        #endif

        bool intersected = false;

        //i hope the compiler doesn't optimize this
        bool tmp = right->intersect(r,is);
        intersected = left->intersect(r,is);



        return tmp | intersected;
    }
    virtual bool intersect_any(Ray& r) const {
        if(!aabb.intersect_any(r))
            return false;
        return right->intersect_any(r) || left->intersect_any(r);
    }
    virtual AABB get_bounds() const {
        return aabb;
    }
};



class BVH: public Primitive {
private:
    int n_partitions = 5;
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::unique_ptr<Primitive> root;

    //template <int axis>
    std::pair<std::vector<int>,std::vector<int>> compute_split(int axis,const std::vector<int>& toSplit, float& cost, const AABB& bounds){

        std::vector<std::vector<int>> partitions(n_partitions);

        for(int i: toSplit){
            bool assigned = false;
            for(int j = 0; j<n_partitions-1; j++){
                if(axis==0){
                if(primitives.at(i)->get_bounds().center().x <
                bounds.min.x + (bounds.max.x-bounds.min.x)*float(j+1)/float(n_partitions)){
                    partitions.at(j).push_back(i);
                    assigned = true;
                }
                } else if(axis ==1) {
                if(primitives.at(i)->get_bounds().center().y <
                bounds.min.y + (bounds.max.y-bounds.min.y)*float(j+1)/float(n_partitions)){
                    partitions.at(j).push_back(i);
                    assigned = true;
                }
                } else {
                if(primitives.at(i)->get_bounds().center().z <
                bounds.min.z + (bounds.max.z-bounds.min.z)*float(j+1)/float(n_partitions)){
                    partitions.at(j).push_back(i);
                    assigned = true;
                }
                }
                if(assigned){
                    break;
                }

            }
            if(!assigned){
                partitions.at(n_partitions-1).push_back(i);
            }
        }
        std::vector<float> partition_cost(n_partitions-1);
        std::vector<AABB> partition_bounds(n_partitions);

        for(int i = 0; i<n_partitions; i++){
            partition_bounds.at(i) = AABB();
            for(int j = 0; j<partitions.at(i).size(); j++){
                partition_bounds.at(i).unite(primitives.at(partitions.at(i).at(j))->get_bounds());
            }
        }
        for(int i = 0; i<n_partitions-1; i++){

            AABB left;
            AABB right;
            int rc = 0;
            int lc = 0;
            for(int j = 0; j<n_partitions; j++){
                for(int k = 0; k<partitions.at(j).size(); k++){
                    if(j<=i){
                        lc+=1;
                        left.unite(primitives.at(partitions.at(j).at(k))->get_bounds());
                    }
                    else {
                        rc+=1;
                        right.unite(primitives.at(partitions.at(j).at(k))->get_bounds());
                    }
                }

            }

            partition_cost.at(i) =  ((lc==0 || rc== 0)?
            std::numeric_limits<float>::infinity() :left.area()
            + right.area());
        }

        int min = 0;
        float min_cost = partition_cost.at(0);
        for(int i = 1; i<n_partitions-1; i++){
            if(partition_cost.at(i)<min_cost){
                min_cost = partition_cost.at(i);
                min = i;
            }
        }
        cost = min_cost;



        std::pair<std::vector<int>,std::vector<int>> ret;
        for(int i = 0; i<=min; i++){
            ret.first.insert(ret.first.end(), partitions.at(i).begin(), partitions.at(i).end());
        }
        for(int i = min+1; i<partitions.size(); i++){
            ret.second.insert(ret.second.end(),partitions.at(i).begin(), partitions.at(i).end());
        }


        return ret;

    }

    std::unique_ptr<Primitive> build_node(std::vector<int>& prims){
        if(prims.size()<4){
            auto g = std::make_unique<SimpleGroup>();
            for(int i = 0; i<prims.size(); i++){
                g->add(primitives.at(prims.at(i)));
            }
            return g;
        }

        AABB bounds;
        AABB bounds_centroid;

        //TODO: compute bounds faster
        for(int i = 0; i<prims.size(); i++){
            AABB tmp = primitives.at(prims.at(i))->get_bounds();
            bounds.unite(tmp);
            bounds_centroid.unite(tmp.center());
        }

        std::vector<int> left;
        std::vector<int> right;
        float c1 = 0.0;
        float c2 = 0.0;
        float c3 = 0.0;
        auto p1 = compute_split(0,prims, c1,bounds_centroid);
        auto p2 = compute_split(1,prims, c2,bounds_centroid);
        auto p3 = compute_split(2,prims, c3,bounds_centroid);

        if(c1<c2){
            if(c1<c3){
                left = p1.first;
                right = p1.second;
            } else if(c2<c3){
                left = p2.first;
                right = p2.second;
            } else {
                left = p3.first;
                right = p3.second;
            }
        } else if(c2<c3){
            left = p2.first;
            right = p2.second;
        } else {

            left = p3.first;
            right = p3.second;
        }



        return std::make_unique<BVHNode>(build_node(left),build_node(right),bounds);
    }

public:
    void add(std::shared_ptr<Primitive> primitive){
        primitives.push_back(primitive);
    }


    void build(){
        auto prims = std::vector<int>(primitives.size());
        for(int i = 0; i<primitives.size(); i++){
            prims.at(i) = i;
        }
        root = build_node(prims);

    }

    bool intersect(Ray& r, Intersection* is) const {
        return root->intersect(r,is);
    }
    bool intersect_any(Ray& r) const {
        return root->intersect_any(r);
    }

    AABB get_bounds() const {
        return root->get_bounds();
    }
};


}

#endif