#include "Triangle.h"
#include <algorithm>
#include <cmath>

Triangle::Triangle()
{
    shared_material = nullptr;
}

float Triangle::intersect(Ray r)
{
    // Extract the vertices of a triangle
    Cartesian3 A = verts[0].Point();
    Cartesian3 B = verts[1].Point();
    Cartesian3 C = verts[2].Point();

    // Calculate the normal vector of a triangle
    Cartesian3 normal = (B - A).cross(C - A);

    // Calculate the intersection point of the ray with the triangular plane
    float t = (A - r.origin).dot(normal) / r.direction.dot(normal);

    // If t is less than zero, it means that the intersection point of the ray and the triangle plane is behind the ray origin.
    if (t < 0.0f) {
        return -1.0f;
    }

    // Calculate intersection point P
    Cartesian3 intersectionPoint = r.origin + t * r.direction;

    // Check if intersection point is inside triangle using half plane algorithm
    Cartesian3 edge1 = B - A;
    Cartesian3 edge2 = C - B;
    Cartesian3 edge3 = A - C;

    Cartesian3 normal1 = edge1.cross(intersectionPoint - A);
    Cartesian3 normal2 = edge2.cross(intersectionPoint - B);
    Cartesian3 normal3 = edge3.cross(intersectionPoint - C);

    // Make sure the intersection point is on the correct triangular half-plane
    if (normal1.dot(normal) >= 0.0f && normal2.dot(normal) >= 0.0f && normal3.dot(normal) >= 0.0f) {
        return t; // Within the triangle, returns the ray parameter t
    }

    return -1.0f; // not within triangle
}


Cartesian3 Triangle::Baricentric(Cartesian3 o){

    // Get the three vertices of the triangle
    Cartesian3 A = verts[0].Point();
    Cartesian3 B = verts[1].Point();
    Cartesian3 C = verts[2].Point();

    // Calculate the normal vector and area of a triangle
    Cartesian3 normal = (B - A).cross(C - A);
    float areaABC = normal.length() * 0.5f;

    // Calculate the vectors from the intersection point to the vertices of two triangles and calculate their respective normal vectors
    Cartesian3 normalA = (o - B).cross(o - C);
    float areaA = normalA.length() * 0.5f;

    Cartesian3 normalB = (o - A).cross(o - C);
    float areaB = normalB.length() * 0.5f;

    Cartesian3 normalC = (o - A).cross(o - B);
    float areaC = normalC.length() * 0.5f;

    // Calculate center of gravity coordinates
    float alpha = areaA / areaABC;
    float beta = areaB / areaABC;
    float gamma = areaC / areaABC;

    // Returns the coordinates
    return Cartesian3(alpha, beta, gamma);

}

Homogeneous4 Triangle::CalculateBlinnPhong(Homogeneous4 lightPosition, Homogeneous4 lightColor, Cartesian3 bc, bool shadow)
{
    // Get material properties from MTL file
    Homogeneous4 ambientColor = shared_material->ambient;  // ambient light color
    Homogeneous4 emissiveColor = shared_material->emissive;  // Self-illuminating color

    // If in shadow, return result early
    if (shadow)
        return ambientColor + emissiveColor;

    // Get material properties from MTL file
    Homogeneous4 diffuseColor = shared_material->diffuse;  // Diffuse light color
    Homogeneous4 specularColor = shared_material->specular;  // Specular light color
    float shininess = shared_material->shininess;  // Gloss

    // Normal vectors at three vertices (m_Normals contains the normal vector at each vertex)
    Cartesian3 normalA = normals[0].Vector();
    Cartesian3 normalB = normals[1].Vector();
    Cartesian3 normalC = normals[2].Vector();

    // Calculate the interpolation normal vector using barycentric coordinate interpolation
    Cartesian3 interpolatedNormal = (bc.x * normalA + bc.y * normalB + bc.z * normalC).unit();

    // Calculate light direction
    Cartesian3 pointOnTriangle = bc.x * verts[0].Point() + bc.y * verts[1].Point() + bc.z * verts[2].Point();
    Cartesian3 lightDirection = (lightPosition.Point() - pointOnTriangle).unit();

    // Calculate diffuse intensity
    float diffuseIntensity = std::max(0.0f, interpolatedNormal.dot(lightDirection));

    // Calculate specular reflection intensity
    Cartesian3 viewDirection = (Cartesian3(0, 0, 1) - pointOnTriangle).unit();
    Cartesian3 halfVector = (lightDirection + viewDirection).unit();
    float specularIntensity = std::pow(std::max(0.0f, interpolatedNormal.dot(halfVector)), shininess);

    // Calculate final color
   /* Homogeneous4 resultColor = lightColor
                                   .modulate(diffuseIntensity * diffuseColor)
                                   .modulate(specularIntensity * specularColor) + ambientColor + emissiveColor;*/

    // Quadratic Attenuation
    float Kc = 1.0f, Kl = 0.09f, Kq = 0.4f;
    float distanceFromLight = (lightPosition.Point() - pointOnTriangle).length();
    float attenuation = float(1.0)/(Kc + Kl * distanceFromLight + Kq * (distanceFromLight * distanceFromLight));



    Homogeneous4 resultColor = lightColor.modulate(diffuseIntensity * diffuseColor) * attenuation +
                                   lightColor.modulate(specularIntensity * specularColor) * attenuation +
                                    ambientColor * attenuation + emissiveColor;

    //resultColor =  ambientColor;
    return resultColor;
}

