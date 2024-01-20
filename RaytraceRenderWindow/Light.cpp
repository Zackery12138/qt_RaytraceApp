//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
////////////////////////////////////////////////////////////////////////

#include "Light.h"
#include <math.h>

Light::Light(LightType type,Homogeneous4 color,Homogeneous4 pos, Homogeneous4 dir, Homogeneous4 tan1, Homogeneous4 tan2)
{
    this->type = type;
    this->lightPosition = pos;
    this->lightDirection = dir;
    this->tangent1 = tan1;
    this->tangent2 = tan2;
    this->lightColor = color;
    enabled = false;
}
Homogeneous4 Light::GetPositionCenter()
{
    return lightPosition;
}

Homogeneous4 Light::SampleAreaLight(){
    float u = -0.5f + (float(rand()) / float(RAND_MAX));
    float v = -0.5f + (float(rand()) / float(RAND_MAX));
    Homogeneous4 pointOnPlane = lightPosition + u * tangent1 + v * tangent2;
    pointOnPlane.w = 1.0f;
    return pointOnPlane;
}

