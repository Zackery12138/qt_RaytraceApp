#include "scene.h"


Scene::Scene(std::vector<ThreeDModel> *texobjs,RenderParameters *renderp)
{
    objects = texobjs;
    rp = renderp;
    Cartesian3 ambient = Cartesian3(0.5f,0.5f,0.5f);
    Cartesian3 diffuse = Cartesian3(0.5f,0.5f,0.5f);
    Cartesian3 specular = Cartesian3(0.5f,0.5f,0.5f);
    Cartesian3 emissive = Cartesian3(0,0,0);
    float shininess = 1.0f;
    default_mat = new Material(ambient,diffuse,specular,emissive,shininess);
}

void Scene::updateScene(){
    triangles.clear(); //clear the list so it can populated again
    for(int i = 0; i < int(objects->size()); i++){
        typedef  unsigned int uint;
        ThreeDModel obj = objects->at(uint(i));
        for(uint face = 0; face < obj.faceVertices.size(); face++){
            for(uint triangle = 0; triangle < obj.faceVertices[face].size() - 2; triangle++){
                Triangle t;
                for (uint vertex = 0; vertex < 3; vertex++) {
                    uint faceVertex = 0;
                    if (vertex != 0)
                        faceVertex = triangle + vertex;
                    //this is our vertex before any transformations. (world space)
                    Homogeneous4 v = Homogeneous4(obj.vertices[obj.faceVertices[face][faceVertex]].x,
                                             obj.vertices[obj.faceVertices[face][faceVertex]].y,
                                             obj.vertices[obj.faceVertices[face][faceVertex]].z);

                    //v = getModelView() * v;
                    t.verts[vertex] = v;

                    Homogeneous4 n = Homogeneous4(obj.normals[obj.faceNormals[face][faceVertex]].x,
                                             obj.normals[obj.faceNormals[face][faceVertex]].y,
                                             obj.normals[obj.faceNormals[face][faceVertex]].z,
                                             0.0f);

                    // n = getModelView() * n;
                    t.normals[vertex] = n;
                    Cartesian3 tex = Cartesian3(obj.textureCoords[obj.faceTexCoords[face][faceVertex]].x,
                                           obj.textureCoords[obj.faceTexCoords[face][faceVertex]].y,
                                           0.0f);
                    t.uvs[vertex] = tex;
                    t.colors[vertex] = Cartesian3( 0.7f, 0.7f, 0.7f);
               }
               if(obj.material== nullptr) {
                    t.shared_material = default_mat;
               }
               else {
                    t.shared_material = obj.material;
               }
               triangles.push_back(t);
            }
        }
    }
}
