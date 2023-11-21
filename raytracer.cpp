/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Leyu Xu
 * *************************
*/
#include "raytracer.h"
#include "Renderable.h"
#include "Sphere.cpp"
#include "Triangle.cpp"
// #ifdef WIN32
//   #include <windows.h>
// #endif

// #if defined(WIN32) || defined(linux)
//   #include <GL/gl.h>
//   #include <GL/glut.h>
// #elif defined(__APPLE__)
//   #include <OpenGL/gl.h>
//   #include <GLUT/glut.h>
// #endif

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #ifdef WIN32
//   #define strcasecmp _stricmp
// #endif

// #include "glm/glm.hpp"

// // #include <imageIO.h>
// #include "png++/png.hpp"
// #include <math.h>
// #include <vector>
// #include <random>

// #define PI 3.14159265

// #define MAX_TRIANGLES 20000
// #define MAX_SPHERES 100
// #define MAX_LIGHTS 100

char * filename = NULL;

// //different display modes
// #define MODE_DISPLAY 1
// #define MODE_JPEG 2

int mode = MODE_DISPLAY;

// //you may want to make these smaller for debugging purposes
// #define WIDTH 640
// #define HEIGHT 480

// //the field of view of the camera
// #define fov 60

unsigned char buffer[HEIGHT][WIDTH][3];



// struct Vertex
// {
//   double position[3];
//   double color_diffuse[3];
//   double color_specular[3];
//   double normal[3];
//   double shininess;
// };

// struct Triangle
// {
//   Vertex v[3];
// };

// struct Sphere
// {
//   double position[3];
//   double color_diffuse[3];
//   double color_specular[3];
//   double shininess;
//   double radius;
// };

struct Light
{
  double position[3];
  double color[3];
};

// struct IntersectData
// {
//     glm::vec3 intersectPoint;
//     glm::vec3 intersectNormal;

//     glm::vec3 interpolatedDiffuseColor;
//     glm::vec3 interpolatedSpecularColor;
//     float interpolatedShininess = 0.0f;
//     float t = FLT_MAX;
// };
Renderable* renderables[MAX_SPHERES + MAX_TRIANGLES];
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
std::vector<Light> subdividedLights;
double ambient_light[3];

// array of screen rays normalized directions
std::vector<glm::vec3> normalizedScreenRaysDirections;

// record pixel colors calculated
unsigned char allPixels[WIDTH][HEIGHT][3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// const float eps = 0.0005f;
// const float eps_accurate = 0.0000001f;

// ****************************** PROGRAM PARAMETERS *************************************
const unsigned int SSAA_Coefficient = 1;

const bool UseSoftShadow = false;
const unsigned int numOfLightSubdivisions = 100;
const float areaLightRadius = 3.0f;

const float recursiveReflection_lambda = 0.3f; // indicates how much the reflection contributes to the lighting
const int recursiveDepth = 0;
// ***************************************************************************************


unsigned char superScaledAllPixels[WIDTH * SSAA_Coefficient][HEIGHT * SSAA_Coefficient][3];

// helper function to convert degrees to radians
float degreesToRadians(float degrees)
{
    return degrees * (PI / 180.0f);
}

// bool isNearlyZero(float val)
// {
//     return abs(val) < eps;
// }

bool areNearlyEqual(float a, float b)
{
    return abs(a - b) < eps;
}

// function to do soft shadow
void SubdivideLightSources()
{
    if (UseSoftShadow)
    {
        std::default_random_engine gen;
        std::uniform_real_distribution<float> distribution(-areaLightRadius, areaLightRadius);

        for (int i = 0; i < num_lights; ++i)
        {
            for (unsigned int j = 0; j < numOfLightSubdivisions; ++j)
            {
                Light randomPoint;

                float randomOffset = distribution(gen);

                glm::vec3 randDirection(distribution(gen), distribution(gen), distribution(gen));
                glm::normalize(randDirection);

                randomPoint.position[0] = lights[i].position[0] + randDirection.x * randomOffset;
                randomPoint.position[1] = lights[i].position[1] + randDirection.y * randomOffset;
                randomPoint.position[2] = lights[i].position[2] + randDirection.z * randomOffset;

                float intensityR = lights[i].color[0] / (float)numOfLightSubdivisions;
                float intensityG = lights[i].color[1] / (float)numOfLightSubdivisions;
                float intensityB = lights[i].color[2] / (float)numOfLightSubdivisions;

                randomPoint.color[0] = intensityR;
                randomPoint.color[1] = intensityG;
                randomPoint.color[2] = intensityB;


                subdividedLights.push_back(randomPoint);
            }
        }
    }
    else
    {
        for (int i = 0; i < num_lights; ++i)
        {
            subdividedLights.push_back(lights[i]);
        }
    }
}

// function to generate all rays shooting to the screen
void generateAllRaysFromCOP()
{
    float aspectRatio = (float)WIDTH / HEIGHT;

    // calculate 3 corners of the screen
    glm::vec3 topLeftScreenCorner;
    topLeftScreenCorner.x = -aspectRatio * tanf(degreesToRadians(fov / 2.0f));
    topLeftScreenCorner.y = tanf(degreesToRadians(fov / 2.0f));
    topLeftScreenCorner.z = -1.0f;

    glm::vec3 topRightScreenCorner;
    topRightScreenCorner.x = -topLeftScreenCorner.x;
    topRightScreenCorner.y = topLeftScreenCorner.y;
    topRightScreenCorner.z = -1.0f;

    glm::vec3 bottomRightScreenCorner;
    bottomRightScreenCorner.x = -topLeftScreenCorner.x;
    bottomRightScreenCorner.y = -topLeftScreenCorner.y;
    bottomRightScreenCorner.z = -1.0f;

    // offset for each pixel
    glm::vec3 pixelWidthOffset = (topRightScreenCorner - topLeftScreenCorner) / ((float)WIDTH * SSAA_Coefficient);
    glm::vec3 pixelHeightOffset = (bottomRightScreenCorner - topRightScreenCorner) / ((float)HEIGHT * SSAA_Coefficient);

    for (int i = 0; i < WIDTH * SSAA_Coefficient; ++i)
    {
        for (int j = HEIGHT * SSAA_Coefficient - 1; j >= 0; --j)
        {
            glm::vec3 ray = topLeftScreenCorner + (float)i * pixelWidthOffset + (float)j * pixelHeightOffset + pixelWidthOffset / 2.0f + pixelHeightOffset / 2.0f;
            ray = glm::normalize(ray);
            normalizedScreenRaysDirections.push_back(ray);
        }
    }
}

// function to test ray-sphere intersection, returns true if intersect and provides relevant data
// bool rayIntersectSphere(const glm::vec3& ray_o, const glm::vec3& ray_d, const Sphere& sphere, IntersectData& data)
// {
//     // if the radius of sphere is zero, no intersection
//     if (sphere.radius < eps)
//     {
//         return false;
//     }

//     float b = 2.0f * (ray_d.x * (ray_o.x - sphere.position[0]) + ray_d.y * (ray_o.y - sphere.position[1]) + ray_d.z * (ray_o.z - sphere.position[2]));
//     float c = powf(ray_o.x - sphere.position[0], 2.0f) + powf(ray_o.y - sphere.position[1], 2.0f) + powf(ray_o.z - sphere.position[2], 2.0f) - powf(sphere.radius, 2.0f);

//     float delta = powf(b, 2.0f) - 4.0f * c;

//     if (delta < 0.0f)
//     {
//         return false;
//     }

//     float t0 = (-b + sqrtf(delta)) / 2.0f;
//     float t1 = (-b - sqrtf(delta)) / 2.0f;

//     // if both t are less or equal to zero, no intersection
//     if (t0 < eps && t1 < eps)
//     {
//         return false;
//     }

//     float t_intersect = 0.0f;

//     if (t0 < eps)
//     {
//         t_intersect = t1;
//     }
//     else if (t1 < eps)
//     {
//         t_intersect = t0;
//     }
//     else
//     {
//         t_intersect = std::min(t0, t1);
//     }

//     // fill out intersection data
//     data.intersectPoint = ray_o + ray_d * t_intersect;
//     glm::vec3 sphereCenter(sphere.position[0], sphere.position[1], sphere.position[2]);
//     data.intersectNormal = (data.intersectPoint - sphereCenter) / (float)sphere.radius;

//     // negate if ray originates inside the sphere
//     if (glm::length(ray_o - sphereCenter) < (float)sphere.radius)
//     {
//         data.intersectNormal = -data.intersectNormal;
//     }

//     data.t = t_intersect;

//     data.interpolatedDiffuseColor = glm::vec3(sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]);
//     data.interpolatedSpecularColor = glm::vec3(sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);
//     data.interpolatedShininess = (float)sphere.shininess;

//     return true;
// }

// function to test ray triangle intersection, returns true if intersect and provide relevant data
// bool rayIntersectTriangle(const glm::vec3& ray_o, const glm::vec3& ray_d, const Triangle& triangle, IntersectData& data)
// {
//     // test if ray and plane are parallel
//     glm::vec3 v1(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
//     glm::vec3 v2(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
//     glm::vec3 v3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

//     glm::vec3 planeNormal = glm::normalize(glm::cross(v2 - v1, v3 - v1));

//     // if ray is parallel to the plane
//     float nDotd = glm::dot(planeNormal, ray_d);
//     if (abs(nDotd) < eps)
//     {
//         return false;
//     }

//     // calculate ray plane intersection
//     float planeCoefficient_d = -(glm::dot(planeNormal, v1));
//     float t = -(glm::dot(planeNormal, ray_o) + planeCoefficient_d) / nDotd;

//     // intersection behind ray origin
//     if (t <= eps)
//     {
//         return false;
//     }

//     // calculate plane intersection point
//     glm::vec3 I = ray_o + t * ray_d;

//     //try to 2D projection onto different planes and calculate barycentric coordinates
//     float areaV1V2V3 = -1.0f;
//     float areaV1V2I = -1.0f;
//     float areaV2V3I = -1.0f;
//     float areaV1IV3 = -1.0f;

//     float weightV1V2I = -1.0f;
//     float weightV2V3I = -1.0f;
//     float weightV1IV3 = -1.0f;

//     // try xy plane
//     if (!isNearlyZero(glm::dot(planeNormal, glm::vec3(0.0f, 0.0f, 1.0f))))
//     {
//         areaV1V2V3 = 0.5f * ((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
//         areaV1V2I = 0.5f * ((v2.x - v1.x) * (I.y - v1.y) - (I.x - v1.x) * (v2.y - v1.y));
//         areaV2V3I = 0.5f * ((v2.x - I.x) * (v3.y - I.y) - (v3.x - I.x) * (v2.y - I.y));
//         areaV1IV3 = 0.5f * ((I.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (I.y - v1.y));
//     }
//     // try xz plane
//     else if (!isNearlyZero(glm::dot(planeNormal, glm::vec3(0.0f, 1.0f, 0.0f))))
//     {
//         areaV1V2V3 = 0.5f * ((v2.x - v1.x) * (v3.z - v1.z) - (v3.x - v1.x) * (v2.z - v1.z));
//         areaV1V2I = 0.5f * ((v2.x - v1.x) * (I.z - v1.z) - (I.x - v1.x) * (v2.z - v1.z));
//         areaV2V3I = 0.5f * ((v2.x - I.x) * (v3.z - I.z) - (v3.x - I.x) * (v2.z - I.z));
//         areaV1IV3 = 0.5f * ((I.x - v1.x) * (v3.z - v1.z) - (v3.x - v1.x) * (I.z - v1.z));
//     }
//     // has to be yz plane
//     else
//     {
//         areaV1V2V3 = 0.5f * ((v2.y - v1.y) * (v3.z - v1.z) - (v3.y - v1.y) * (v2.z - v1.z));
//         areaV1V2I = 0.5f * ((v2.y - v1.y) * (I.z - v1.z) - (I.y - v1.y) * (v2.z - v1.z));
//         areaV2V3I = 0.5f * ((v2.y - I.y) * (v3.z - I.z) - (v3.y - I.y) * (v2.z - I.z));
//         areaV1IV3 = 0.5f * ((I.y - v1.y) * (v3.z - v1.z) - (v3.y - v1.y) * (I.z - v1.z));
//     }

//     // test if the intersection is in triangle, reject if we get negative weight
//     weightV1V2I = areaV1V2I / areaV1V2V3;
//     if (weightV1V2I < 0.0f)
//     {
//         return false;
//     }

//     weightV2V3I = areaV2V3I / areaV1V2V3;
//     if (weightV2V3I < 0.0f)
//     {
//         return false;
//     }

//     weightV1IV3 = areaV1IV3 / areaV1V2V3;
//     if (weightV1IV3 < 0.0f)
//     {
//         return false;
//     }

//     // fill out intersection data
//     glm::vec3 v1_normal(triangle.v[0].normal[0], triangle.v[0].normal[1], triangle.v[0].normal[2]);
//     glm::vec3 v2_normal(triangle.v[1].normal[0], triangle.v[1].normal[1], triangle.v[1].normal[2]);
//     glm::vec3 v3_normal(triangle.v[2].normal[0], triangle.v[2].normal[1], triangle.v[2].normal[2]);

//     glm::vec3 v1_diffuse(triangle.v[0].color_diffuse[0], triangle.v[0].color_diffuse[1], triangle.v[0].color_diffuse[2]);
//     glm::vec3 v2_diffuse(triangle.v[1].color_diffuse[0], triangle.v[1].color_diffuse[1], triangle.v[1].color_diffuse[2]);
//     glm::vec3 v3_diffuse(triangle.v[2].color_diffuse[0], triangle.v[2].color_diffuse[1], triangle.v[2].color_diffuse[2]);

//     glm::vec3 v1_specular(triangle.v[0].color_specular[0], triangle.v[0].color_specular[1], triangle.v[0].color_specular[2]);
//     glm::vec3 v2_specular(triangle.v[1].color_specular[0], triangle.v[1].color_specular[1], triangle.v[1].color_specular[2]);
//     glm::vec3 v3_specular(triangle.v[2].color_specular[0], triangle.v[2].color_specular[1], triangle.v[2].color_specular[2]);

//     float v1_shininess = triangle.v[0].shininess;
//     float v2_shininess = triangle.v[1].shininess;
//     float v3_shininess = triangle.v[2].shininess;

//     // interpolate
//     data.intersectPoint = I;
//     data.intersectNormal = glm::normalize(v1_normal * weightV2V3I + v2_normal * weightV1IV3 + v3_normal * weightV1V2I);
//     data.interpolatedDiffuseColor = v1_diffuse * weightV2V3I + v2_diffuse * weightV1IV3 + v3_diffuse * weightV1V2I;
//     data.interpolatedSpecularColor = v1_specular * weightV2V3I + v2_specular * weightV1IV3 + v3_specular * weightV1V2I;
//     data.interpolatedShininess = v1_shininess * weightV2V3I + v2_shininess * weightV1IV3 + v3_shininess * weightV1V2I;
//     data.t = t;

//     return true;
// }

// function to test segment-sphere intersection, returns true if the segment intersects the sphere
bool segmentIntersectSphere(const glm::vec3& a, const glm::vec3& b, const Sphere& sphere)
{
    IntersectData data;
    glm::vec3 ab = b - a;
    glm::vec3 ray_d = glm::normalize(ab);

    if (sphere.intersectRay(a, ray_d, data)) //rayIntersectSphere(a, ray_d, sphere, data)
    {
        glm::vec3 aToIntersectionPoint = data.t * ray_d;

        float lengthI = glm::length(aToIntersectionPoint);
        float lengthAB = glm::length(ab);

        // if the intersection point is not in the segment
        if (lengthI >= lengthAB)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        return false;
    }
}

bool segmentIntersectTriangle(const glm::vec3& a, const glm::vec3& b, const Triangle& triangle)
{
    IntersectData data;
    glm::vec3 ab = b - a;
    glm::vec3 ray_d = glm::normalize(ab);

    if (triangle.intersectRay(a, ray_d, data))//(rayIntersectTriangle(a, ray_d, triangle, data))
    {
        glm::vec3 aToIntersectionPoint = data.t * ray_d;

        float lengthI = glm::length(aToIntersectionPoint);
        float lengthAB = glm::length(ab);

        // if the intersection point is not in the segment
        if (lengthI >= lengthAB)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    else
    {
        return false;
    }
}

// function to do recursive reflection
glm::vec3 recursiveRayTrace(const glm::vec3& ray_o, const glm::vec3& ray_d, int depth)
{
    // base case
    if (depth <= 0)
    {
        return glm::vec3(0.0f, 0.0f, 0.0f);
    }
    else
    {
        // go through all spheres in scene and find the intersection with the smallest t
        IntersectData data;
        IntersectData tempData;
        for (int k = 0; k < num_spheres; ++k)
        {
            if (spheres[k].intersectRay(ray_o, ray_d, tempData)) //(rayIntersectSphere(ray_o, ray_d, spheres[k], tempData))
            {
                // if we see a smaller t update intersection data
                if (tempData.t < data.t)
                {
                    data = tempData;
                }
            }
        }

        // go through all triangles in scene and find the intersection with the smallest t
        IntersectData data_tri;
        IntersectData tempData_tri;
        for (int l = 0; l < num_triangles; ++l)
        {
            if (triangles[l].intersectRay(ray_o, ray_d, tempData_tri))//(rayIntersectTriangle(ray_o, ray_d, triangles[l], tempData_tri))
            {
                if (tempData_tri.t < data_tri.t)
                {
                    data_tri = tempData_tri;
                }
            }
        }

        // choose between smallest t sphere and triangle
        if (data_tri.t < data.t)
        {
            data = data_tri;
        }

        // if we find intersection
        if (data.t < FLT_MAX)
        {
            glm::vec3 pixelColor(0.0f, 0.0f, 0.0f);

            // add ambient color once
            pixelColor.x += ambient_light[0];
            pixelColor.y += ambient_light[1];
            pixelColor.z += ambient_light[2];

            // go through each light source and cast a shadow ray
            for (unsigned int m = 0; m < subdividedLights.size(); ++m)
            {
                glm::vec3 a = data.intersectPoint;
                glm::vec3 b(subdividedLights[m].position[0], subdividedLights[m].position[1], subdividedLights[m].position[2]);

                bool isInShadow = false;

                // test whether we're shadowed by a sphere
                for (int n = 0; n < num_spheres; ++n)
                {
                    if (segmentIntersectSphere(a, b, spheres[n]))
                    {
                        isInShadow = true;
                        break;
                    }
                }

                // test whether we're shadowed by a triangle
                for (int g = 0; g < num_triangles; ++g)
                {
                    if (segmentIntersectTriangle(a, b, triangles[g]))
                    {
                        isInShadow = true;
                        break;
                    }
                }

                // if we're not in shadow in terms of this light source, do phong shading
                if (!isInShadow)
                {
                    glm::vec3 lightPosition(subdividedLights[m].position[0], subdividedLights[m].position[1], subdividedLights[m].position[2]);
                    glm::vec3 L = glm::normalize(lightPosition - data.intersectPoint);
                    float LdotN = glm::dot(L, data.intersectNormal);

                    // clamp dot product
                    if (LdotN < 0.0f)
                    {
                        LdotN = 0.0f;
                    }

                    glm::vec3 R = glm::normalize(-glm::reflect(L, data.intersectNormal));
                    glm::vec3 V = glm::normalize(glm::vec3(0.0f, 0.0f, 0.0f) - data.intersectPoint);
                    float RdotV = glm::dot(R, V);

                    // clamp dot product
                    if (RdotV < 0.0f)
                    {
                        RdotV = 0.0f;
                    }

                    // calculate color for each channel due to this light source, and add it to the final color

                    // red
                    float r =
                        subdividedLights[m].color[0] * (data.interpolatedDiffuseColor.x * LdotN + data.interpolatedSpecularColor.x * powf(RdotV, data.interpolatedShininess));

                    pixelColor.x += r;
                    // clamp if necessary
                    if (pixelColor.x > 1.0f)
                    {
                        pixelColor.x = 1.0f;
                    }

                    // green
                    float g =
                        subdividedLights[m].color[1] * (data.interpolatedDiffuseColor.y * LdotN + data.interpolatedSpecularColor.y * powf(RdotV, data.interpolatedShininess));

                    pixelColor.y += g;
                    // clamp if necessary
                    if (pixelColor.y > 1.0f)
                    {
                        pixelColor.y = 1.0f;
                    }

                    // blue
                    float b =
                        subdividedLights[m].color[2] * (data.interpolatedDiffuseColor.z * LdotN + data.interpolatedSpecularColor.z * powf(RdotV, data.interpolatedShininess));

                    pixelColor.z += b;
                    // clamp if necessary
                    if (pixelColor.z > 1.0f)
                    {
                        pixelColor.z = 1.0f;
                    }
                }
            }

            glm::vec3 newDirection = glm::normalize(glm::reflect(ray_d, data.intersectNormal));
            glm::vec3 reflectedColor = recursiveRayTrace(data.intersectPoint, newDirection, depth - 1);
            glm::vec3 finalColor = (1.0f - recursiveReflection_lambda) * pixelColor + recursiveReflection_lambda * reflectedColor;

            return finalColor;
        }
        else
        {
            return glm::vec3(0.0f, 0.0f, 0.0f);
        }
    }
}

// function to process the entire scene
void processScene()
{
    SubdivideLightSources();

    int rayIndex = 0;

    // go through each pixel
    for (int i = 0; i < WIDTH * SSAA_Coefficient; ++i)
    {
        for (int j = 0; j < HEIGHT * SSAA_Coefficient; ++j)
        {
            // go through all spheres in scene and find the intersection with the smallest t
            IntersectData data;
            IntersectData tempData;
            for (int k = 0; k < num_spheres; ++k)
            {
                // if (rayIntersectSphere(glm::vec3(0.0f, 0.0f, 0.0f), normalizedScreenRaysDirections[rayIndex], spheres[k], tempData))
                if(spheres[k].intersectRay(glm::vec3(0.0f, 0.0f, 0.0f), normalizedScreenRaysDirections[rayIndex], tempData))
                {
                    // if we see a smaller t update intersection data
                    if (tempData.t < data.t)
                    {
                        data = tempData;
                    }
                }
            }

            // go through all triangles in scene and find the intersection with the smallest t
            IntersectData data_tri;
            IntersectData tempData_tri;
            for (int l = 0; l < num_triangles; ++l)
            {
                //if (rayIntersectTriangle(glm::vec3(0.0f, 0.0f, 0.0f), normalizedScreenRaysDirections[rayIndex], triangles[l], tempData_tri))
                if(triangles[l].intersectRay(glm::vec3(0.0f, 0.0f, 0.0f), normalizedScreenRaysDirections[rayIndex],tempData_tri))
                {
                    if (tempData_tri.t < data_tri.t)
                    {
                        data_tri = tempData_tri;
                    }
                }
            }

            // choose between smallest t sphere and triangle
            if (data_tri.t < data.t)
            {
                data = data_tri;
            }


            // if we find intersection
            if (data.t < FLT_MAX)
            {
                glm::vec3 pixelColor(0.0f, 0.0f, 0.0f);

                // add ambient color once
                pixelColor.x += ambient_light[0];
                pixelColor.y += ambient_light[1];
                pixelColor.z += ambient_light[2];

                // go through each light source and cast a shadow ray
                for (unsigned int m = 0; m < subdividedLights.size(); ++m)
                {
                    glm::vec3 a = data.intersectPoint;
                    glm::vec3 b(subdividedLights[m].position[0], subdividedLights[m].position[1], subdividedLights[m].position[2]);

                    bool isInShadow = false;

                    // test whether we're shadowed by a sphere
                    for (int n = 0; n < num_spheres; ++n)
                    {
                        if (segmentIntersectSphere(a, b, spheres[n]))
                        {
                            isInShadow = true;
                            break;
                        }
                    }

                    // test whether we're shadowed by a triangle
                    for (int g = 0; g < num_triangles; ++g)
                    {
                        if (segmentIntersectTriangle(a, b, triangles[g]))
                        {
                            isInShadow = true;
                            break;
                        }
                    }

                    // if we're not in shadow in terms of this light source, do phong shading
                    if (!isInShadow)
                    {
                        glm::vec3 lightPosition(subdividedLights[m].position[0], subdividedLights[m].position[1], subdividedLights[m].position[2]);
                        glm::vec3 L = glm::normalize(lightPosition - data.intersectPoint);
                        float LdotN = glm::dot(L, data.intersectNormal);

                        // clamp dot product
                        if (LdotN < 0.0f)
                        {
                            LdotN = 0.0f;
                        }

                        glm::vec3 R = glm::normalize(-glm::reflect(L, data.intersectNormal));
                        glm::vec3 V = glm::normalize(glm::vec3(0.0f, 0.0f, 0.0f) - data.intersectPoint);
                        float RdotV = glm::dot(R, V);

                        // clamp dot product
                        if (RdotV < 0.0f)
                        {
                            RdotV = 0.0f;
                        }

                        // calculate color for each channel due to this light source, and add it to the final color
                        
                        // red
                        float r =
                            subdividedLights[m].color[0] * (data.interpolatedDiffuseColor.x * LdotN + data.interpolatedSpecularColor.x * powf(RdotV, data.interpolatedShininess));

                        pixelColor.x += r;
                        // clamp if necessary
                        if (pixelColor.x > 1.0f)
                        {
                            pixelColor.x = 1.0f;
                        }

                        // green
                        float g = 
                            subdividedLights[m].color[1] * (data.interpolatedDiffuseColor.y * LdotN + data.interpolatedSpecularColor.y * powf(RdotV, data.interpolatedShininess));

                        pixelColor.y += g;
                        // clamp if necessary
                        if (pixelColor.y > 1.0f)
                        {
                            pixelColor.y = 1.0f;
                        }

                        // blue
                        float b = 
                            subdividedLights[m].color[2] * (data.interpolatedDiffuseColor.z * LdotN + data.interpolatedSpecularColor.z * powf(RdotV, data.interpolatedShininess));

                        pixelColor.z += b;
                        // clamp if necessary
                        if (pixelColor.z > 1.0f)
                        {
                            pixelColor.z = 1.0f;
                        }
                    }
                }

                if (recursiveDepth <= 0)
                {
                    superScaledAllPixels[i][j][0] = (unsigned char)(255.0f * pixelColor.x);
                    superScaledAllPixels[i][j][1] = (unsigned char)(255.0f * pixelColor.y);
                    superScaledAllPixels[i][j][2] = (unsigned char)(255.0f * pixelColor.z);
                }
                else
                {
                    glm::vec3 newDirection = glm::normalize(glm::reflect(normalizedScreenRaysDirections[rayIndex], data.intersectNormal));
                    glm::vec3 reflectedColor = recursiveRayTrace(data.intersectPoint, newDirection, recursiveDepth);
                    glm::vec3 finalColor = (1.0f - recursiveReflection_lambda) * pixelColor + recursiveReflection_lambda * reflectedColor;

                    superScaledAllPixels[i][j][0] = (unsigned char)(255.0f * finalColor.x);
                    superScaledAllPixels[i][j][1] = (unsigned char)(255.0f * finalColor.y);
                    superScaledAllPixels[i][j][2] = (unsigned char)(255.0f * finalColor.z);
                }
                
            }
            // else assign the background color
            else
            {
                superScaledAllPixels[i][j][0] = 255;
                superScaledAllPixels[i][j][1] = 255;
                superScaledAllPixels[i][j][2] = 255;
            }

            ++rayIndex;
        }
    }

    // scale down for anti-aliasing
    for (int i = 0, m = 0; i < WIDTH && m < WIDTH * SSAA_Coefficient - 1; ++i, m += SSAA_Coefficient)
    {
        for (int j = 0, n = 0; j < HEIGHT && n < HEIGHT * SSAA_Coefficient - 1; ++j, n += SSAA_Coefficient)
        {
            float divider = (float)SSAA_Coefficient * (float)SSAA_Coefficient;
            
            float r = 0.0f;
            float g = 0.0f;
            float b = 0.0f;

            for (int u = 0; u < SSAA_Coefficient; ++u)
            {
                for (int v = 0; v < SSAA_Coefficient; ++v)
                {
                    r += (float)superScaledAllPixels[m + u][n + v][0];
                    g += (float)superScaledAllPixels[m + u][n + v][1];
                    b += (float)superScaledAllPixels[m + u][n + v][2];
                }
            }

            allPixels[i][j][0] = (unsigned char)(r / divider);

            allPixels[i][j][1] = (unsigned char)(g / divider);

            allPixels[i][j][2] = (unsigned char)(b / divider);
        }
    }
}

//MODIFY THIS FUNCTION
void draw_scene()
{
    // generate all rays from COP
    generateAllRaysFromCOP();

    // process the scene
    processScene();

  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      plot_pixel(x, y, allPixels[x][y][0], allPixels[x][y][1], allPixels[x][y][2]);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

// void save_jpg()
// {
//   printf("Saving JPEG file: %s\n", filename);

//   ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
//   if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
//     printf("Error in Saving\n");
//   else 
//     printf("File saved Successfully\n");
// }

void save_png()
{
    png::image< png::rgb_pixel > image(WIDTH, HEIGHT);
    for (png::uint_32 y = 0; y < image.get_height(); ++y)
    {
        for (png::uint_32 x = 0; x < image.get_width(); ++x)
        {
            image[y][x] = png::rgb_pixel(buffer[x][y][0], buffer[x][y][1], buffer[x][y][2]);
            // non-checking equivalent of image.set_pixel(x, y, ...);
        }
    }
    image.write("rgb.png");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
//   Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");
      Sphere* newSphere= new Sphere;
      parse_doubles(file,"pos:",newSphere->position);
      parse_rad(file,&newSphere->radius);
      parse_doubles(file,"dif:",newSphere->color_diffuse);
      parse_doubles(file,"spe:",newSphere->color_specular);
      parse_shi(file,&newSphere->shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = *newSphere;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_png();
    int random=0;
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

