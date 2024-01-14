#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <float.h>
#include "ThreeDModel.h"
#include "RenderParameters.h"
#include "Triangle.h"
#include "Ray.h"


class Scene
{
public:
    struct CollisionInfo
    {
        Triangle tri;
        float t;
        Cartesian3 CollisionPoint;
        bool Valid(){return t > float(1e-6) && t < FLT_MAX;}
    };

    std::vector<ThreeDModel>* objects;
    RenderParameters* rp;
    std::vector<Triangle> triangles;
    Material *default_mat;


public:
    Scene(std::vector<ThreeDModel> *texobjs, RenderParameters *renderp);
    void updateScene();
    Matrix4 getModelView();

    CollisionInfo closestTriangle(Ray r);

};

#endif // SCENE_H
