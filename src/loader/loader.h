#ifndef CPPPT_LOADER
#define CPPPT_LOADER

#include <renderer/renderer.h>
#include <scene.h>

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <texture/image_texture.h>
#include <image/rgb_image.h>
#include <bxdf/emissive_material.h>

namespace cpppt
{
namespace Loader{


    /*private*/
    namespace {
        struct LoaderData {
            Scene* s;
            SceneData* sd;
            RenderSettings* rs;
            std::unordered_map<std::string,int> loaded_materials;

        };
        void parse_float(const std::string& s, float* v){
            *v = atof((s.substr(6, s.length()-7)).c_str());
        }


        void parse_rgb(const std::string& s, float* r, float* g, float* b){
            int c1 = s.find(",");

            std::string s2 = s.substr(4,c1-4);
            *r = atof(s2.c_str());
            s2 = s.substr(c1+1);
            int c2 = s2.find(",");
            *g = atof((s2.substr(0,c2)).c_str());
            *b = atof((s2.substr(c2+1,s2.length()-c2-2)).c_str());
        }

        template <bool bw, bool srgb>
        std::shared_ptr<Texture> get_texture(std::string directory, std::string tex){
            if(tex.substr(0,3) == "rgb"){
                float r,g,b;
                parse_rgb(tex,&r,&g,&b);
                if(!bw){
                    Vec3 v(r,g,b);
                    /*if(srgb){
                        v = pixel_ops::srgb2linear(v);
                    }*/
                    return std::make_shared<ConstantTexture<Vec3>>(v);
                } else {
                    /*if(srgb){
                        r = pixel_ops::srgb2linear(r);
                    }*/
                    return std::make_shared<ConstantTexture<float>>(r);
                }
            } else if(tex.substr(0,5) == "float"){
                float r;
                std::stringstream ss(tex);
                parse_float(tex,&r);
                if(!bw){
                    Vec3 v(r);
                    if(srgb){
                        v = pixel_ops::srgb2linear(v);
                    }
                    return std::make_shared<ConstantTexture<Vec3>>(v);
                } else {
                    if(srgb){
                        r = pixel_ops::srgb2linear(r);
                    }
                    return std::make_shared<ConstantTexture<float>>(r);
                }
            } else if(tex.substr(0,5) == "image"){
                tex = tex.substr(6,tex.length()-7);
                auto color_space = ImageDef::LINEAR;
                if(srgb){
                    color_space = ImageDef::SRGB;
                }
                if(!bw){
                    return std::make_shared<ImageTexture<Vec3>>(
                        std::make_shared<Image<Vec3>>(directory+tex,color_space));
                }else{
                    return std::make_shared<ImageTexture<float>>(
                        std::make_shared<Image<float>>(directory+tex,color_space));
                }


            } else {
                throw std::runtime_error("Unknown type of tex: " + tex);
                return nullptr;
            }
        }

        std::string trim(
            const std::string& s,
            const std::string& whitespace = " \n\t\r")
        {
            int i = s.find_first_not_of(whitespace);
            int j = s.find_last_not_of(whitespace);
            return s.substr(i,j-i+1);
        }

        std::pair<std::string,std::string> split_filename (const std::string& str)
        {
            std::pair<std::string,std::string> ret;
            int idx =str.find_last_of("/\\");
            ret.first = str.substr(0,idx) + "/";
            ret.second = str.substr(idx+1);
            return ret;
        }

        inline void parse_camera(LoaderData* ld, std::istringstream& iss){

            float tan_fovy;
            Vec3 origin;
            Vec3 at;
            Vec3 up;
            iss>>origin.x>>origin.y>>origin.z>>at.x>>at.y>>at.z>>up.x>>up.y>>up.z>>tan_fovy;
            ld->s->camera = std::make_shared<CameraPerspective>(
                Vec2i(ld->rs->resX,ld->rs->resY), tan_fovy, origin, normalized(at-origin), up
            );
        }
        inline void parse_light(LoaderData* ld, std::istringstream& iss){
            std::string type;
            iss>>type;
            if(type=="point"){
                Vec3 pos;
                std::string tex;
                float intensity;
                iss>>pos.x>>pos.y>>pos.z>>tex>>intensity;
                //TODO: directory
                ld->s->light->add(std::make_shared<PointLight>(pos,get_texture<false,true>("",tex),intensity));
            } else {
                throw std::runtime_error("Unknown light type: "+ type);
                return;
            }
        }


        inline void parse_mat(LoaderData* ld, std::string dir, std::string filename){
            auto it = ld->loaded_materials.find(filename);
            if(it!=ld->loaded_materials.end()){
                ld->sd->materials.push_back(ld->sd->materials.at(it->second));
            } else {
                std::ifstream iss(dir+filename);
                std::string word;
                iss>>word;
                if(word == "DISNEY"){
                    iss>>word;
                    std::shared_ptr<Texture> albedo = get_texture<false, true>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> spec = get_texture<true, true>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> roughness = get_texture<true, false>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> metal = get_texture<true, false>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> transparent = get_texture<true,true>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> ior = get_texture<true, true>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> normal = get_texture<false, false>(dir,word);
                    iss>>word;
                    std::shared_ptr<Texture> alpha = get_texture<true, false>(dir,word);
                    iss>>word;

                    std::shared_ptr<StandardMaterial> mat = std::make_shared<StandardMaterial>(
                        albedo,normal,spec,roughness, metal,ior,transparent,alpha
                    );

                    ld->loaded_materials.insert(std::pair<std::string,int>
                    (filename, static_cast<int>( ld->sd->materials.size()))
                    );
                    ld->sd->materials.push_back(mat);
                } else if(word=="EMISSION") {

                    float intensity = 1.0;
                    iss>>word;
                    std::shared_ptr<Texture> color = get_texture<false, true>(dir,word);
                    /*TODO: check*//*
                    Vec3 col = color->sample(Vec3(0.8,0.2,0.0));
                    std::cout<<col.x<<" "<<col.y<<" "<<col.z<<std::endl;
                    */
                    iss>>word;

                    intensity = atof((word.substr(6,word.length()-7)).c_str());

                    std::shared_ptr<EmissiveMaterial> mat = std::make_shared<EmissiveMaterial>(
                        color, intensity
                    );
                    ld->loaded_materials.insert(std::pair<std::string,int>
                        (filename, static_cast<int>(ld->sd->materials.size()))
                    );
                    ld->sd->materials.push_back(mat);

                } else {
                    throw std::runtime_error("Wrong bsdf type: "+ word);
                }
            }

        }

        inline void parse_obj(LoaderData* ld, std::string dir, std::string filename){
            int n_v;
            int n_f;
            std::string type;
            std::ifstream ifs(dir+filename);

            std::string line;
            if(!std::getline(ifs,line)){
                throw std::runtime_error("Insuficient number of lines when reading file");
            }
            std::istringstream iss(line);
            iss>>n_v>>n_f>>type;
            if(type!= "TRIANGLE_MESH"){
                throw std::runtime_error("Unknown mesh type: "+ type);
            }
            ld->sd->meshes.push_back(std::make_shared<Mesh>(n_v,n_f));

            for(int i = 0; i<n_v; i++){
                if(!std::getline(ifs,line)){
                    throw std::runtime_error("Insuficient number of lines when reading verts");
                }

                std::istringstream iss(line);
                Vec3 p;
                Vec2 uv;
                Vec3 normal;
                Vec3 tangent;
                iss>>p.x>>p.y>>p.z>>uv.x>>uv.y>>normal.x>>normal.y>>normal.z>>tangent.x>>tangent.y>>tangent.z;
                Vec3 bitangent = cross(normal,tangent);

                ld->sd->meshes.at(ld->sd->meshes.size()-1)->set_vertex_data(i,p,normal,uv,bitangent);
            }

            for(int i = 0; i<n_f; i++){
                if(!std::getline(ifs,line)){
                    throw std::runtime_error("Insuficient number of lines when reading faces");
                }
                std::istringstream iss(line);
                int a,b,c;
                iss>>a>>b>>c;
                ld->sd->meshes.at(ld->sd->meshes.size()-1)->set_face_data(i,a,b,c);
            }
        }

    }


    inline void load_render_settings(RenderSettings* rs, std::string filename){
        /*First line: renderer name*/
        /*second line: resx resy */
        /*third line: n samples*/

        std::ifstream ifs(filename);

        std::string line;
        while(std::getline(ifs,line)){
            line = line.substr(0, line.find("#"));
            std::istringstream iss(line);
            std::string type;
            iss>>type;
            if(type == "renderer"){
                std::string rtype;
                iss>>rtype;
                if(rtype=="debug"){
                    rs->renderer = RenderSettings::RendererType::DEBUG;
                } else if(rtype=="lighttracer"){
                    rs->renderer = RenderSettings::RendererType::LIGHTTRACER;
                } else if(rtype=="pathtracer"){
                    rs->renderer = RenderSettings::RendererType::PATHTRACER;
                } else {
                    throw std::runtime_error("Unknown renderer type: " + rtype);
                }
            } else if (type == "spp") {
                iss>>rs->spp;
            } else if (type == "resolution") {
                iss>>rs->resX>>rs->resY;
            } else if (type!="") {
                throw std::runtime_error("Unknown render setting: " + type);
            }

        }



    }


    inline void load_scene(
        Scene* s,
        SceneData* sd,
        RenderSettings* rs,
        std::string filename){


        s->primitive = std::make_shared<BVH>();
        s->light = std::make_shared<LightGroup>();
        auto dir_file = split_filename(filename);
        LoaderData ld;
        ld.s = s;
        ld.sd = sd;
        ld.rs = rs;
        std::ifstream ifs(filename);
        std::string line;
        while(std::getline(ifs,line)){
            line = line.substr(0, line.find("#"));
            line = trim(line);

            std::istringstream iss(line);
            std::string type;
            iss>>type;
            if(type == "camera"){
                parse_camera(&ld,iss);
            } else if (type == "obj") {
                std::string obj;
                std::string mat;
                iss>>obj>>mat;
                parse_obj(&ld, dir_file.first, obj);
                parse_mat(&ld, dir_file.first, mat);

            } else if (type == "light") {
                parse_light(&ld,iss);
            } else if (type!="") {
                throw std::runtime_error("Unknown scene obj: " + type);
            }
        }

        for(int i = 0; i<sd->meshes.size(); i++){
            for(int j = 0; j<sd->meshes.at(i)->n_triangles(); j++){
                auto helper = std::make_shared<PrimitiveLeaf>(
                    std::make_shared<Triangle>(sd->meshes.at(i).get(),j),
                    sd->materials.at(i)
                );

                s->primitive->add(helper);

                if(sd->materials.at(i)->is_emissive()){
                    s->light->add(std::make_shared<ShapeLight>(helper));
                }

            }
        }
        std::cout<<"Done loading..."<<std::endl;

        std::cout<<"building accel..."<<std::endl;

        s->primitive->build();
        //s->light->build();
    }

}


} // namespace cppt


#endif