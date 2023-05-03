import bpy
from mathutils import Vector
from math import tan
from math import pi
import os
import shutil

# note: we assume:
# every mesh has identity transform
# only triangles
# everything has a uv map
# each mesh has only one material
# modifiers applied
    
path = '~/cpp-pt/scenes/ring/'
filename = "scene"

def save_tex(node):
    print("saving_tex:",node)
    try:
        image_path = bpy.path.abspath(node.links[0].from_node.image.filepath)
        new_image_path = os.path.join(path,"images",os.path.basename(image_path))
        shutil.copy(image_path, new_image_path)
        return "image("+os.path.join("images",os.path.basename(image_path))+")"
    except Exception as e:
        print(e)
        v = node.default_value
        if(type(v) == type(0.5)):
            return "float("+str(v)+")"
        else:
            return "rgb("+str(v[0])+","+str(v[1])+","+str(v[2])+")"

def save_normal(node):
    print("saving_normal:",node)
    try:
        image_path = bpy.path.abspath(node.links[0].from_node.inputs["Color"].links[0].from_node.image.filepath)
        new_image_path = os.path.join(path,"images",os.path.basename(image_path))
        shutil.copy(image_path, new_image_path)
        return "image("+os.path.join("images",os.path.basename(image_path))+")"
    except Exception as e:
        print(e)
        return "rgb(0.5,0.5,1.0)"

def save_env(f):
    node = bpy.data.worlds[0].node_tree.nodes["World Output"].inputs["Surface"].links[0].from_node
    tex = save_tex(node.inputs["Color"])
    if(not tex.startswith("image")):
        return
    else:
        print("environment", tex, "2048 1024",file=f)

def save_material(material, filename):
    
    print("saving_mat:",material)
    node = material.node_tree.nodes["Material Output"].inputs['Surface'].links[0].from_node;
    albedo = "rgb(0.8,0.0,0.8)"
    spec = "float(0.5)"
    roughness = "float(0.5)"
    metal = "float(0.0)"
    transparent = "float(0.0)"
    ior = "float(1.3)"
    normal = "rgb(0.5,0.5,1.0)"
    alpha = "float(1.0)"
    
    
    if(node.name == 'Principled BSDF'):
        albedo = save_tex(node.inputs["Base Color"])
        spec = save_tex(node.inputs["Specular"])
        roughness = save_tex(node.inputs["Roughness"])
        metal = save_tex(node.inputs["Metallic"])
        transparent = save_tex(node.inputs["Transmission"])
        ior = save_tex(node.inputs["IOR"])
        normal = save_normal(node.inputs["Normal"])
        alpha = save_tex(node.inputs["Alpha"])
    elif(node.name == "Emission"):
        f = open(os.path.join(path,"materials",filename),"w")
        color = save_tex(node.inputs["Color"])
        strength = save_tex(node.inputs["Strength"])
        print("EMISSION",file=f)
        print(color,file =f )
        print(strength,file =f )
        f.close()
        return
    
    
    f = open(os.path.join(path,"materials",filename),"w")
    print("DISNEY",file=f)
    print(albedo,file =f )
    print(spec,file =f )
    print(roughness,file =f )
    print(metal,file =f )
    print(transparent,file =f )
    print(ior,file =f )
    print(normal,file =f )
    print(alpha,file =f )
    f.close()

    #order: albedo spec roughness metal transparent ior normal alpha

def save_object(obj,fdes):
    print("saving_objs:",obj.data.name)
    material_name = obj.data.materials[0].name
    
    print("obj", "models/" + obj.name + ".ptobj", "materials/" + material_name+".ptmat", file = fdes)
    save_material(obj.data.materials[0], material_name + ".ptmat")
    
    obj.data.calc_loop_triangles()
    obj.data.calc_tangents()
    f2 = open(os.path.join(path,"models",obj.name+".ptobj"), "w")
    uv_layer_data = obj.data.uv_layers[0].data;
    
    seen_vert = [[]]*len(obj.data.vertices)
    
    verts = []
    uvs = []
    normals = []
    faces = []
    tangents = []
    #verts = [None]*len(uv_layer_data.uv)
    #uvs = uv_layer_data.uv
    #normals = [None]*len(uv_layer_data.uv)

    #print(len(obj.data.polygons), len(uv_layer_data),file = f2)
    helper_i = 0
    for tri in obj.data.polygons:
        face = []
        for id in tri.loop_indices:
            v_id = obj.data.loops[id].vertex_index;


            
            helper_i +=1
            if helper_i %100 == 0:
                print(helper_i)
            #print(v_id)
            #if vertex has been seen_verten, just write the corresponding uv idx
            exists_vert = False;
            #if len(seen_vert[v_id]) > 0:
            #    for v in seen_vert[v_id]:
            #        if uvs[v] == uv_layer_data[id].uv:
            #            face.append(v)
            #            exists_vert = True
            #            break;
            if not exists_vert:
                new_v_id = len(verts)
                face.append(new_v_id)
                verts.append(obj.data.vertices[v_id].co)
                normals.append(obj.data.vertices[v_id].normal)
                tangents.append(obj.data.loops[id].tangent)
                #tood: bitan sign?
                uvs.append(uv_layer_data[id].uv)
                seen_vert[v_id].append(new_v_id)

        faces.append(face)
    
    print(len(verts), len(faces),"TRIANGLE_MESH",file=f2)
    
    for i in range(len(verts)):
        print(str(verts[i].x), str(verts[i].y),str(verts[i].z),
            str(uvs[i].x), str(uvs[i].y),
            str(normals[i].x), str(normals[i].y),str(normals[i].z),
            str(tangents[i].x), str(tangents[i].y),str(tangents[i].z),
            file=f2)
            
    for i in range(len(faces)):
        print(str(faces[i][0]),str(faces[i][1]),str(faces[i][2]), file = f2)
            
            

if not os.path.exists(path):
    os.makedirs(path)
    
if not os.path.exists(os.path.join(path,"images")):
    os.makedirs(os.path.join(path,"images"))

if not os.path.exists(os.path.join(path,"models")):
    os.makedirs(os.path.join(path,"models"))


if not os.path.exists(os.path.join(path,"materials")):
    os.makedirs(os.path.join(path,"materials"))

def save_light(ob,f):
    print("light ", end="",file=f)
    if(ob.data.type == 'POINT'):
        print("point", ob.location.x,ob.location.y,ob.location.z,
        "rgb("+
        str(ob.data.color.r)+","+
        str(ob.data.color.g)+","+
        str(ob.data.color.b)+")",
        ob.data.energy, file = f)
def save_camera(ob,f):
    loc = ob.location
    at = ob.matrix_world@(Vector([0.0,0.0,-1.0,1.0]))
    up = ob.matrix_world@(Vector([0.0,1.0,0.0,0.0]))
    tan_fovy = ob.data.sensor_width/(ob.data.lens)*0.5
    print("camera", loc.x,loc.y,loc.z,at.x,at.y,at.z,up.x,up.y,up.z,tan_fovy,file = f)


f = open(path+filename+".ptscene","w");
camera_saved = False

save_env(f)

for ob in bpy.data.objects:
    if ob.type == 'MESH':
        save_object(ob,f)
    
    if ob.type == 'LIGHT':
        save_light(ob,f)
    if ob.type == 'CAMERA':
        if not camera_saved:
            save_camera(ob,f)
            camera_saved = True





f.close()