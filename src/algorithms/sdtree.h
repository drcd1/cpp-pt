#ifndef CPPPT_SDTREE
#define CPPPT_SDTREE
#include <math/math.h>
#include <math/sampler.h>

#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <iomanip>

namespace cpppt{

template <typename T>
struct Bounds{
    T min;
    T max;
    T center(){
        return (min+max)*0.5;
    }
    Bounds():min(10e6),max(-10e6){}
    Bounds(const T& a, const T& b):min(a),max(b){}
};




template
<typename Leaf, int d>
class TreeNode{
public:
    TreeNode():isLeaf(true),leaf(){
    }
    TreeNode(Leaf l):isLeaf(true), leaf(l){
    }
    //TreeNode(const TreeNode<Leaf,d>& n) = delete;
   /*
        TreeNode(TreeNode<Leaf,d>&& n):

        children(std::move(n.children)),

    {
    }
*/

    TreeNode& operator=(const TreeNode<Leaf,d>& n){
        isLeaf = n.isLeaf;
        if(isLeaf){
            leaf = Leaf(n.leaf);
            return *this;
        }
        children = new std::array<TreeNode<Leaf, d>,d>;
        for(int i = 0; i<d; i++){
            children->at(i) = n.children->at(i);
        }
        return *this;
    }
    TreeNode(const TreeNode<Leaf,d>& n){
        *this = n;
        /*
        if(isLeaf){
            leaf = Leaf(n.leaf);
            return;
        }
        children = new std::array<TreeNode<Leaf, d>,d>;
        for(int i = 0; i<d; i++){
            children->at(i) = n.children->at(i);
        }
        */
    }

    ~TreeNode(){
        if(!isLeaf){
            delete children;
        }
    }
    bool isLeaf;
    union{
        std::array<TreeNode<Leaf, d>,d>* children;
        Leaf leaf;
    };
};


typedef TreeNode<float,4> DNode;


struct DNodeBounds {
    DNodeBounds(DNode* n, Bounds<Vec2>& b):n(n),b(b){}
    DNode* n;
    Bounds<Vec2> b;
};


struct DTree{
private:

    static float get_radiance_and_prune(DNode* node, float threshold){
        if(node->isLeaf){
            return node->leaf;
        }
        else {
            float rad = 0.0;
            for(int i =0; i<4; i++){
                rad += get_radiance_and_prune(&(node->children->at(i)),threshold);
            }
            if(rad<threshold){
                delete node->children;
                node->isLeaf = true;
                node->leaf = rad;
            }
            return rad;
        }

    }
public:
    //TODO: multiply radiance by ct factor?
    DTree(): root(0.0),count(0),radiance(0){
        leaves.push_back(Leaf({-M_PI,-1},{2.0*M_PI,2}));
    }
    DTree(const DTree& d): root(d.root),count(d.count), radiance(d.radiance), leaves(d.leaves){
       // std::cout<<"Creating dtree: "<<leaves.size()<<std::endl;
       // std::cout<<"From: "<<d.leaves.size()<<std::endl;
    }
    ~DTree(){

    };
    DNode root;

    struct Leaf{
        Leaf(const Vec2& begin, const Vec2& size): begin(begin),size(size){}
        Vec2 begin;
        Vec2 size;
    };
    std::vector<Leaf> leaves = std::vector<Leaf>();

    template <typename DLeafFunction>
    void iterate_leaves(DLeafFunction dlf){
        std::vector<DNodeBounds> stack;
        Bounds<Vec2> b;
        b.min=Vec2(-M_PI,-1.0);
        b.max=Vec2(M_PI,1.0);
        DNodeBounds a(&root,b);
        stack.push_back(DNodeBounds(&root,b));

        while(!stack.empty()){
            DNodeBounds d = stack.back();
            stack.pop_back();
            if(d.n->isLeaf){
                dlf(d);
            } else {
                Bounds<Vec2> b = d.b;
                Vec2 size = (b.max - b.min)*0.5;
                for(int i = 0; i<4; i++){
                    Vec2 begin = Vec2(b.min.x+size.x*(i%2), b.min.y + size.y*(int(i/2)));
                    Bounds<Vec2> b2;
                    b2.min = begin;
                    b2.max = begin+size;
                    if(begin.x<-M_PI-0.01 ||b2.max.x>M_PI+0.01){
                        std::cout<<"somethingwrong"<<std::endl;
                    }
                    stack.push_back(DNodeBounds(&(d.n->children->at(i)), b2));
                }
            }
        }
    }


    //adapt_and_build_leaves

    //iterate_leaves

    void adapt_and_build_leaves(float factor = 1.0){
        if(radiance == 0.0){
            return;
        }
        std::cout<<"Radiance Total: "<<radiance<<std::endl;
        float threshold = radiance*0.01;
        leaves = std::vector<Leaf>();
        std::vector<DNode*> stack;
        std::vector<Leaf> leaf_stack;
        stack.push_back(&root);
        leaf_stack.push_back(Leaf(Vec2(-M_PI, -1),Vec2(2.0*M_PI,2.0)));
        //two pass algorithm: first, prune tree. Then, build leaves and subdivision

        get_radiance_and_prune(&root, threshold);
        float new_rad = 0.0;
        while(!stack.empty()){
            DNode* dn = stack.back();
            stack.pop_back();
            Leaf l = leaf_stack.back();
            leaf_stack.pop_back();
            if(dn->isLeaf){
                float ratio = dn->leaf/threshold;
                if(ratio>1.0){
                //do subdivision
                    //DNode a(std::move(DNode(dn->leaf*0.25)));
                    dn->isLeaf = false;
                    dn->children = new std::array<DNode,4>{
                        DNode(dn->leaf*0.25),
                        DNode(dn->leaf*0.25),
                        DNode(dn->leaf*0.25),
                        DNode(dn->leaf*0.25)
                    };
                } else {
                    new_rad+=dn->leaf;
                    dn->leaf = 0.0;
                    //todo: set radiance to zero?
                    leaves.push_back(l);
                    continue;
                }
            }
            for(int i = 0; i<4;i++){
                stack.push_back(&(dn->children->at(i)));
                leaf_stack.push_back(Leaf(
                    Vec2(l.begin.x+l.size.x*0.5*(i%2), l.begin.y+l.size.y*0.5*(int(i/2))),
                    l.size*0.5
                    )
                );
            }
        }

        std::cout<<"Total radiance: "<<new_rad<<std::endl;


        //TODO: Maybe we can just set everything to zero?
        //this way, the pdf is encoded only in the structure of the tree
        //each leaf node is sampled uniformly
        //and the pdf is 1/(node_area*n_nodes)
        count = 0;
        radiance = 0.0;
        //note: count was never set?
        //radiance = radiance*factor;
    };
    int count;
    float radiance;
};
template
<typename T>
class MovedPtr{
public:
    MovedPtr(){};
    MovedPtr(std::unique_ptr<T> a){
        pt = std::move(a);
    }
    MovedPtr(MovedPtr<T>& a){
        pt = std::move(a.pt);
    }
    MovedPtr operator=(MovedPtr<T>& a){
        pt = std::move(a.pt);
        return this;
    }
    ~MovedPtr() = default;
    std::unique_ptr<T> pt;
};


typedef TreeNode<MovedPtr<DTree>,2> SNode;


    struct SNodeBounds {
        SNodeBounds(SNode* n, Bounds<Vec3>& b, int dim):n(n),b(b),dim(dim){}
        SNode* n;
        Bounds<Vec3> b;
        int dim;
    };


class SDTree{
    bool learn = true;

    Bounds<Vec3> spatial_bounds;
    SNode root;
    int iteration = 1;
    int c = 24000;


    Vec2 get_area(DNode& root, float x, float y) const {
        int depth = 0;

        DNode* node = &root;
        Bounds<Vec2> cb;
        cb.min = Vec2(-M_PI, -1.0);
        cb.max = Vec2(M_PI, 1.0);
        while(true){

            Vec2 mid = cb.center();
            if(node->isLeaf){
                return cb.max-cb.min;
            }
            int child;
            if(x<mid.x){
                if(y<mid.y){
                    child = 0;
                    cb.max = mid;
                } else {
                    child = 2;
                    cb.max.x = mid.x;
                    cb.min.y = mid.y;
                }
            } else {
                if(y<mid.y){
                    child = 1;
                    cb.min.x = mid.x;
                    cb.max.y = mid.y;
                } else {
                    child = 3;
                    cb.min = mid;
                }
            }
            node = &node->children->at(child);

        }

    }

    float& get_leaf(DNode& root, float x, float y) const {
        int depth = 0;

        DNode* node = &root;
        Bounds<Vec2> cb;
        cb.min = Vec2(-M_PI, -1.0);
        cb.max = Vec2(M_PI, 1.0);
        while(true){

            Vec2 mid = cb.center();
            if(node->isLeaf){
                return node->leaf;
            }
            int child;
            if(x<mid.x){
                if(y<mid.y){
                    child = 0;
                    cb.max = mid;
                } else {
                    child = 2;
                    cb.max.x = mid.x;
                    cb.min.y = mid.y;
                }
            } else {
                if(y<mid.y){
                    child = 1;
                    cb.min.x = mid.x;
                    cb.max.y = mid.y;
                } else {
                    child = 3;
                    cb.min = mid;
                }
            }
            node = &node->children->at(child);

        }
    }

    DTree& get_leaf(const Vec3& pt) const {
        int depth = 0;
        const SNode* node = &root;

        Bounds<Vec3> cb = spatial_bounds;

        while(true){
            if(node->isLeaf){
                return *(node->leaf.pt);
            }
            bool right = false;
            float mid;
            switch(depth%3){
                case 0:
                    mid = (cb.min.x + cb.max.x)*0.5;
                    right = pt.x>mid;
                    if(!right){
                        cb.max.x= mid;
                        node = &node->children->at(0);
                    } else {
                        cb.min.x = mid;
                        node = &node->children->at(1);
                    }
                    break;
                case 1:
                    mid = (cb.min.y + cb.max.y)*0.5;
                    right = pt.y>mid;
                    if(!right){
                        cb.max.y= mid;
                        node = &node->children->at(0);
                    } else {
                        cb.min.y = mid;
                        node = &node->children->at(1);
                    }
                    break;
                case 2:
                    mid = (cb.min.z + cb.max.z)*0.5;
                    right = pt.z>mid;
                    if(!right){
                        cb.max.z= mid;
                        node = &node->children->at(0);
                    } else {
                        cb.min.z = mid;
                        node = &node->children->at(1);
                    }
                    break;
                default:
                    throw std::runtime_error("Should not reach");
            }

            depth++;
        }
    }

    float sample_direction(const DTree& leaf, Vec3* wi, Sampler& s) const {
        //TODO: we have two options: either just store the leaves of the tree and do binary search on prob, or
        int sz = leaf.leaves.size();
        int i = (int)(s.sample()*sz);
        float u = s.sample();
        float v = s.sample();
        //this may be very wrong
        DTree::Leaf dtl = leaf.leaves.at(i%sz);
        float area = dtl.size.x*dtl.size.y;

        float phi = dtl.begin.x+u*dtl.size.x;
        float cos_theta = dtl.begin.y + v*dtl.size.y;
        float sin_theta = sqrt(1.0-(cos_theta*cos_theta));
        *wi = Vec3(cos(phi)*sin_theta,sin(phi)*sin_theta,cos_theta);
        //this may be very wrong todo

        float pdf =  1.0/(leaf.leaves.size()*area);
        if(std::isnan(pdf) || pdf==0.0){
            std::cout<<"New Nan!"<<std::endl;
        }
        return pdf;

    }

public:
    void no_learning(){
        learn = false;
    }


    SDTree(const Bounds<Vec3>& spatial_bounds): spatial_bounds(spatial_bounds),
    root(std::make_unique<DTree>())
    {

    }


    DirectionalSample sample(Sampler& s, const Vec3& point) const {
        DTree& leaf = get_leaf(point);
        DirectionalSample ds;
        ds.pdf = sample_direction(leaf, &(ds.wi), s);
        ds.delta = false;
        return ds;
    }

    void update(){
        subdivide(c,iteration);
        iteration++;
    }


    void subdivide(int c,int iteration){
        int threshold = c*int(pow(2,iteration*0.5));

        std::vector<SNode*> stack;
        std::vector<int> subdivided;
        stack.push_back(&root);
        subdivided.push_back(0);
        int depth=0;
        std::cout<<"threshold:"<<threshold<<std::endl;
        while(!stack.empty()){
            SNode* n = stack.back();
            int sub = subdivided.back();
            stack.pop_back();
            subdivided.pop_back();


            if(n->isLeaf){
                if(n->leaf.pt->count>threshold || sub){
                    //how many times to subdivide
                    int to_sub = sub;
                    if(sub==0){
                        //compute how many times we will subdivide:
                        float ratio = float(n->leaf.pt->count)/float(threshold);
                        to_sub = int(log2(ratio))+1;
                        float factor = 1.0/float(to_sub);
                        //we haven't subdivided yet, so we have to build the quadtree
                        // which will be copied through all subdivided
                        // tree
                        n->leaf.pt->adapt_and_build_leaves(factor);
                    }
                    std::array<SNode,2>* a1 = new std::array<SNode,2>{
                        MovedPtr<DTree>(std::make_unique<DTree>(*(n->leaf.pt))),
                        MovedPtr<DTree>(std::make_unique<DTree>(*(n->leaf.pt)))};
                    n->isLeaf = false;
                    n->children = a1;
                    stack.push_back(&(n->children->at(0)));
                    subdivided.push_back(to_sub-1);
                    stack.push_back(&(n->children->at(1)));
                    subdivided.push_back(to_sub-1);
                } else {
                    //we have not subdivided yet, and we won't
                    //but we need to rebuild the quadtree
                    n->leaf.pt->adapt_and_build_leaves(1.0);
                }

            } else {
                stack.push_back(&(n->children->at(0)));
                subdivided.push_back(false);
                stack.push_back(&(n->children->at(1)));
                subdivided.push_back(false);
            }
        }

    }

    void splat(float radiance, Vec3 pos, Vec3 dir){
        DTree& leaf = get_leaf(pos);
        float phi = atan2(dir.y,dir.x);
        float cos_theta = dir.z;
        #pragma omp atomic
        (get_leaf(leaf.root,phi,cos_theta)) += radiance;
        #pragma omp atomic
        leaf.count += 1;
        #pragma omp atomic
        leaf.radiance += radiance;
    }
    float pdf(Vec3 pos, Vec3 dir) const{
        DTree& leaf = get_leaf(pos);
        float phi = atan2(dir.y,dir.x);
        float cos_theta = dir.z;
        Vec2 area = get_area(leaf.root,phi,cos_theta);
        float pdf = (1.0/leaf.leaves.size())*(1.0/(area.x*area.y));

        if(std::isnan(pdf) || pdf == 0.0){
            std::cout<<"more nan"<<std::endl;
        }

        return pdf;

    }
    template <typename SLeafFunction, typename DLeafFunction>
    void iterate_leaves(SLeafFunction slf, DLeafFunction dlf) {



        std::vector<SNodeBounds> stack;



        stack.push_back(SNodeBounds(&root, spatial_bounds,0));

        while(!stack.empty()){
            SNodeBounds d = stack.back();
            stack.pop_back();
            int dim = d.dim;
            if(d.n->isLeaf){
                slf(d);
                d.n->leaf.pt->iterate_leaves(dlf);
            } else {
                Bounds<Vec3> b1 = d.b;
                Bounds<Vec3> b2 = d.b;

                if(dim==0){
                    b1.max.x = (d.b.min.x+d.b.max.x)*0.5;
                    b2.min.x = (d.b.min.x+d.b.max.x)*0.5;
                } else if(dim==1){
                    b1.max.y = (d.b.min.y+d.b.max.y)*0.5;
                    b2.min.y = (d.b.min.y+d.b.max.y)*0.5;
                } else if(dim==2){
                    b1.max.z = (d.b.min.z+d.b.max.z)*0.5;
                    b2.min.z = (d.b.min.z+d.b.max.z)*0.5;
                }

                stack.push_back(SNodeBounds(&(d.n->children->at(0)),b1,(dim+1)%3));
                stack.push_back(SNodeBounds(&(d.n->children->at(1)),b2,(dim+1)%3));
            }
        }
    }

    //notes: DTree is subdivided if node has >1% of DTree rad
    //STree is subdivided if node has > c*(sqrt(2^k)) verts, where c is 12000, k is num iteration

    void debug_log(std::string folder) {
        std::cout<<"opening file helloworld.txt"<<std::endl;
        std::ofstream test("helloworld.txt");
        test<<"Hello World!"<<std::endl;
        test.close();
        std::ofstream ofs("debug_log/"+folder+"/stree.obj");
        std::ofstream ofs2;
        std::ostringstream oss;

        std::vector<Vec3> verts;
        std::vector<std::array<int,4>> faces;
        int n_leaf = 0;
        int n_leaf_leaf = 0;
        int i2 = 0;
        iterate_leaves(
        [&](const SNodeBounds& snb){
            //finish the previous_node:
            if(n_leaf>0){
                ofs2.close();
            }
            n_leaf_leaf = 0;

            std::ostringstream oss2;
            if(n_leaf_leaf>0){
            for(int i = 0; i<n_leaf_leaf; i++){

            }
            ofs2.close();
            }
            oss.str(std::string());
            oss<<std::setfill('0')<<std::setw(10);
            oss<<n_leaf;

            int i = n_leaf;
            oss.str(std::string());
            oss<<std::setfill('0')<<std::setw(10);
            oss<<i;
            ofs<<"o "<<oss.str()<<std::endl;

            ofs<<"v "<<snb.b.min.x<<" "<<snb.b.min.y<<" "<<snb.b.min.z<<std::endl;
            ofs<<"v "<<snb.b.min.x<<" "<<snb.b.min.y<<" "<<snb.b.max.z<<std::endl;
            ofs<<"v "<<snb.b.min.x<<" "<<snb.b.max.y<<" "<<snb.b.min.z<<std::endl;
            ofs<<"v "<<snb.b.min.x<<" "<<snb.b.max.y<<" "<<snb.b.max.z<<std::endl;
            ofs<<"v "<<snb.b.max.x<<" "<<snb.b.min.y<<" "<<snb.b.min.z<<std::endl;
            ofs<<"v "<<snb.b.max.x<<" "<<snb.b.min.y<<" "<<snb.b.max.z<<std::endl;
            ofs<<"v "<<snb.b.max.x<<" "<<snb.b.max.y<<" "<<snb.b.min.z<<std::endl;
            ofs<<"v "<<snb.b.max.x<<" "<<snb.b.max.y<<" "<<snb.b.max.z<<std::endl;


            /*binary code:
            z = 0 : even
            z = 1: odd
            y = 0 : /2 even
            y = 1: /2 odd
            x = 0
            */

            ofs<<"f "<<1+i*8+0<<" "<<1+i*8+2<<" "<<1+i*8+6<<" "<<1+i*8+4<<std::endl;
            ofs<<"f "<<1+i*8+1<<" "<<1+i*8+3<<" "<<1+i*8+7<<" "<<1+i*8+5<<std::endl;
            ofs<<"f "<<1+i*8+0<<" "<<1+i*8+1<<" "<<1+i*8+5<<" "<<1+i*8+4<<std::endl;
            ofs<<"f "<<1+i*8+2<<" "<<1+i*8+3<<" "<<1+i*8+7<<" "<<1+i*8+6<<std::endl;
            ofs<<"f "<<1+i*8+0<<" "<<1+i*8+1<<" "<<1+i*8+3<<" "<<1+i*8+2<<std::endl;
            ofs<<"f "<<1+i*8+4<<" "<<1+i*8+5<<" "<<1+i*8+7<<" "<<1+i*8+6<<std::endl;

            ofs2.open("debug_log/"+ folder+"/"+oss.str() + ".obj");

            n_leaf++;

        },
        [&](const DNodeBounds& snb){
            int i = n_leaf_leaf;

            std::ostringstream oss2;
            oss2<<std::setfill('0')<<std::setw(10)<<"l"<<i<<std::endl;
            ofs2<<"o "<<oss2.str()<<std::endl;
            ofs2<<"v "<<snb.b.min.x<<" "<<snb.b.min.y<<" 0.0\n";
            ofs2<<"v "<<snb.b.max.x<<" "<<snb.b.min.y<<" 0.0\n";
            ofs2<<"v "<<snb.b.max.x<<" "<<snb.b.max.y<<" 0.0\n";
            ofs2<<"v "<<snb.b.min.x<<" "<<snb.b.max.y<<" 0.0\n";
            ofs2<<"f "<<1+i*4<<" "<<1+i*4+1<<" "<<1+i*4+2<<" "<<1+i*4+3<<std::endl;
            n_leaf_leaf++;
        }

        );

        ofs2.close();
        ofs.close();



    }


};


};
#endif