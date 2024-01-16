//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <random>
#include <QTimer>
// include the header file
#include "RaytraceRenderWidget.h"

#define N_THREADS 16
#define N_LOOPS 100
#define N_BOUNCES 5
#define TERMINATION_FACTOR 0.35f

// constructor
RaytraceRenderWidget::RaytraceRenderWidget
        (   
        // the geometric object to show
        std::vector<ThreeDModel>      *newTexturedObject,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // parent widget in visual hierarchy
        QWidget             *parent
        )
    // the : indicates variable instantiation rather than arbitrary code
    // it is considered good style to use it where possible
    : 
    // start by calling inherited constructor with parent widget's pointer
    QOpenGLWidget(parent),
    // then store the pointers that were passed in
    texturedObjects(newTexturedObject),
    renderParameters(newRenderParameters)
    { // constructor
        QTimer *timer = new QTimer(this);
        connect(timer, &QTimer::timeout, this, &RaytraceRenderWidget::forceRepaint);
        timer->start(30);
        scene = new Scene(newTexturedObject, newRenderParameters);
    } // constructor


// destructor
RaytraceRenderWidget::~RaytraceRenderWidget()
    { // destructor
    // empty (for now)
    // all of our pointers are to data owned by another class
    // so we have no responsibility for destruction
    // and OpenGL cleanup is taken care of by Qt
    } // destructor                                                                 

// mouse-handling
void RaytraceRenderWidget::mousePressEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mousePressEvent()
    // store the button for future reference
    int whichButton = int(event->button());
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;

    
    // and we want to force mouse buttons to allow shift-click to be the same as right-click
    unsigned int modifiers = event->modifiers();
    
    // shift-click (any) counts as right click
    if (modifiers & Qt::ShiftModifier)
        whichButton = Qt::RightButton;
    
    // send signal to the controller for detailed processing
    emit BeginScaledDrag(whichButton, x,y);
    } // RaytraceRenderWidget::mousePressEvent()
    
void RaytraceRenderWidget::mouseMoveEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseMoveEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit ContinueScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseMoveEvent()
    
void RaytraceRenderWidget::mouseReleaseEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseReleaseEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit EndScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseReleaseEvent()

// called when OpenGL context is set up
void RaytraceRenderWidget::initializeGL()
    { // RaytraceRenderWidget::initializeGL()
	// this should remain empty
    } // RaytraceRenderWidget::initializeGL()

// called every time the widget is resized
void RaytraceRenderWidget::resizeGL(int w, int h)
    { // RaytraceRenderWidget::resizeGL()
    // resize the render image
    frameBuffer.Resize(w, h);
    } // RaytraceRenderWidget::resizeGL()
    
// called every time the widget needs painting
void RaytraceRenderWidget::paintGL()
    { // RaytraceRenderWidget::paintGL()
    // set background colour to white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    // and display the image
    glDrawPixels(frameBuffer.width, frameBuffer.height, GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer.block);
    } // RaytraceRenderWidget::paintGL()


void RaytraceRenderWidget::Raytrace(){
    raytracingThread = std::thread(&RaytraceRenderWidget::RaytraceThread, this);
    raytracingThread.detach();
}

void RaytraceRenderWidget::RaytraceThread(){
    scene->updateScene();
    frameBuffer.clear(RGBAValue(0.0f, 0.0f, 0.0f, 1.0f));

    #pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < frameBuffer.height; j++){
        for(int i = 0; i < frameBuffer.width; i++){
            Ray ray = calculateRay(i, j, !renderParameters->orthoProjection);

            //Homogeneous4 color(i/float(frameBuffer.height), j/float(frameBuffer.width),0);
            auto color = traceAndShadeWithRay(ray, N_BOUNCES,1.0f);




            //Gamma correction...
            float gamma = 2.2f;
            //We already calculate everything in float, so we just do gamma correction
            //before putting it interger format;
            color.x = pow(color.x, 1/gamma);
            color.y = pow(color.y, 1/gamma);
            color.z = pow(color.z, 1/gamma);
            frameBuffer[j][i] = RGBAValue(color.x * 255.0f, color.y * 255.0f, color.z * 255.0f, 255.0f);
        }
    }
    std::cout << "Done!" << std::endl;
}

void RaytraceRenderWidget::forceRepaint(){
    update();
}

Ray RaytraceRenderWidget::calculateRay(int pixelx, int pixely, bool perspective){
    // Convert pixel coordinates to normalized device coordinates (NDC)
    float ndcX = (float(2.0) * pixelx / frameBuffer.width - float(1.0));
    float ndcY = (float(2.0) * pixely / frameBuffer.height - float(1.0));

    // Calculate aspect ratio
    float aspect = float(1.0) * frameBuffer.width / frameBuffer.height;

    float viewX, viewY;

    // Adjust aspect ratio
    if (aspect > 1)
    {
        viewX = ndcX * aspect;
        viewY = ndcY;
    }
    else
    {
        viewX = ndcX;
        viewY = ndcY / aspect;
    }

    // Define rays based on orthogonal or perspective projection
    // perspective projection
    if (perspective)
    {
        // Map 2D coordinates to points on the image plane
        Cartesian3 positionOnImage(viewX, viewY, 0);

        // Defines the starting point of the ray at (0, 0, 1), directed towards the image plane
        Cartesian3 rayOrigin(0, 0, 1);
        Cartesian3 direction = (positionOnImage - rayOrigin).unit();

        Ray perspectiveRay(rayOrigin, direction);
        return perspectiveRay;
    }
    // Orthographic projection
    else
    {
        Cartesian3 rayOrigin(viewX, viewY, 0);
        Cartesian3 direction(0, 0, -1);
        Ray orthoRay(rayOrigin, direction);
        return orthoRay;
    }
}

Homogeneous4 RaytraceRenderWidget::traceAndShadeWithRay(Ray ray, int depth, float curIor)
{
    Homogeneous4 color(0, 0, 0, 0);

    if (depth == 0)
    {
        return color;
    }

    // Try to find triangles in the scene where rays collide
    auto collisionInfo = scene->closestTriangle(ray);

    // A collision occurs!!!
    if (collisionInfo.Valid())
    {
        auto collisionPoint = collisionInfo.CollisionPoint;
        auto triangle = collisionInfo.tri;
        auto bc = triangle.Baricentric(collisionPoint);
        auto normOut = (bc.x * triangle.normals[0] + bc.y * triangle.normals[1] + bc.z * triangle.normals[2]).Vector().unit();

        // Recursive ray tracing
        // Get material properties
        auto reflectivity = triangle.shared_material->reflectivity;
        auto surfaceIOR = triangle.shared_material->indexOfRefraction;
        auto transparency = triangle.shared_material->transparency;

        // interpolation rendering
        if (renderParameters->interpolationRendering)
        {
            color = Homogeneous4(abs(normOut.x), abs(normOut.y), abs(normOut.z), 1);
        }
        // Blinn-Phong lighting model
        else if (renderParameters->phongEnabled)
        {

            color = Homogeneous4(0.0, 0.0, 0.0, 1.0);
            for (const auto& light : renderParameters->lights)
            {
                // Apply model view matrix
                auto lightPosition = scene->getModelView() * light->GetPositionCenter();

                // Calculate shadow effectiveness
                bool shadowValid = renderParameters->shadowsEnabled && calculateShadowValidity(triangle, collisionPoint, lightPosition);

                // Calculate Blinn-Phong lighting
                color = color + triangle.CalculateBlinnPhong(lightPosition, light->GetColor(), bc, shadowValid);
            }
            //reflection
            if(reflectivity > float(0.0) && renderParameters->reflectionEnabled){
                Ray mirrorRay = reflectRay(ray, normOut, collisionPoint);


                return reflectivity * traceAndShadeWithRay(mirrorRay, depth-1, curIor) + (1 - reflectivity) * color;
            }
            //refraction
            if(transparency > float(0.0) && renderParameters->refractionEnabled){
                Ray refractedRay = refractRay(ray, normOut, collisionPoint, curIor, surfaceIOR);

                return traceAndShadeWithRay(refractedRay, depth - 1, surfaceIOR);
            }

            // Fresnel effect
            if(reflectivity > float(0.0) && transparency > float(0.0) && renderParameters-> fresnelRendering){
                float cosTheta = ray.direction.unit().dot(normOut);

                float rTheta = fresnel_schilick(cosTheta, curIor, surfaceIOR);
                float updatedReflectivity = reflectivity * rTheta;
                float updatedReTransparency = transparency * (1 - rTheta);
                Ray mirrorRay = reflectRay(ray, normOut, collisionPoint);
                color = color + updatedReflectivity * traceAndShadeWithRay(mirrorRay, depth - 1, curIor);
                Ray refractedRay = refractRay(ray, normOut, collisionPoint, curIor, surfaceIOR);
                color = color + updatedReTransparency * traceAndShadeWithRay(refractedRay, depth - 1, surfaceIOR);
            }


        }
        else
        {
            // interpolated color
            color = bc.x * triangle.colors[0] + bc.y * triangle.colors[1] + bc.z * triangle.colors[2];
        }
    }

    return color;
}

bool RaytraceRenderWidget::calculateShadowValidity(const Triangle& triangle, const Cartesian3& collisionPoint, const Homogeneous4& lightPosition)
{
    // Compute shadow rays
    auto collisionPointNormal = triangle.normals[0].Vector();//Get the normals of a triangle
    // avoid shadow acne
    auto shadowRayOrigin = collisionPoint + 1e-4 * collisionPointNormal;// Calculate the starting point of the shadow ray, slightly offset from the collision point along the normal direction
    auto shadowRayDirection = (lightPosition.Point() - shadowRayOrigin).unit();// Calculate the direction of the shadow ray, pointing to the light source position
    Ray shadowRay(shadowRayOrigin, shadowRayDirection);// Create shadow ray

    //Try to find the object in the scene that the light collides with
    Scene::CollisionInfo shadowCollisionInfo = scene->closestTriangle(shadowRay);

    if (shadowCollisionInfo.Valid())
    {
        // Check if the shadow is valid.
        auto shadow2Collision = shadowCollisionInfo.CollisionPoint - collisionPoint;
        auto collision2Light = lightPosition.Point() - shadowCollisionInfo.CollisionPoint;

        // Shadows have no effect if the collision point is a light source
        if (collision2Light.length() < 1e-6 || shadow2Collision.length() < 1e-6)
        {
            return false;
        }

        // If the vector from shadowRayOrigin to shadowHitPoint is equal to
        // The vectors from shadowHitPoint to lightPosition have the same direction
        // then it is a valid shadow
        return shadow2Collision.dot(collision2Light) > 0;
    }

    return true;
}


Ray RaytraceRenderWidget::reflectRay(Ray &ray, Cartesian3 &normOut, Cartesian3 collisionPoint){
    //calculate the reflected ray's direction
    auto reflectDirection = (ray.direction - 2 * (ray.direction.dot(normOut)) * normOut).unit();
    //construct the reflected ray
    Ray reflectedRay(collisionPoint, reflectDirection);
    return reflectedRay;

}
Ray RaytraceRenderWidget::refractRay(Ray &ray, Cartesian3 normOut, Cartesian3 collisionPoint, float inIOR, float outIOR){
    float cosi = ray.direction.unit().dot(normOut);
    float etai = inIOR, etat = outIOR;
    Cartesian3 n = normOut;
    Cartesian3 refractedDirection;
    if(cosi < 0){
        cosi = -cosi;
    }else {
        // Light exits from the inside, swapping the refractive index, flipping the normal.
        std::swap(etai, etat);
        n = Cartesian3(-normOut.x, -normOut.y, -normOut.z);
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    if(k < 0){
        //total internal reflection
        return reflectRay(ray, n, collisionPoint);
        //std::cout << "total internal reflection , we do nothing cuttently"<<std::endl;
    }
    refractedDirection = (eta * ray.direction.unit() + (eta * cosi - std::sqrt(k)) * n).unit();

    Ray refractedRay(collisionPoint, refractedDirection);
    return refractedRay;


}



float RaytraceRenderWidget::fresnel_schilick(float cosTheta, float ior1, float ior2){
    float r0 = (ior1 - ior2) / (ior1 + ior2);
    r0 = r0 * r0;
    float rTheta = r0 + (float(1.0) - r0) * float(pow((float(1.0) - cosTheta), 5));
    return rTheta;
}

