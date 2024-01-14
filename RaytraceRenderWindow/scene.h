#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "ThreeDModel.h"
#include "RenderParameters.h"
#include "triangle.h"
class Scene
{
public:
    std::vector<ThreeDModel>* objects;
    RenderParameters* rp;
    std::vector<Triangle> triangles;
    Material *default_mat;


    Scene(std::vector<ThreeDModel> *texobjs, RenderParameters *renderp);
    void updateScene();

};

#endif // SCENE_H
