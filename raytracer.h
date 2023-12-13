#ifndef RAYTRACER_HEADER
#define RAYTRACER_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include "libs/glm/glm.hpp"

#include "libs/lodepng/lodepng.h"
#include <math.h>
#include <vector>
#include <random>
#include <iostream>

#define PI 3.14159265

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

// char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

// int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60

// unsigned char buffer[HEIGHT][WIDTH][3];

const float eps = 0.0005f;
const float eps_accurate = 0.0000001f;

struct IntersectData
{
    public:
    glm::vec3 intersectPoint;
    glm::vec3 intersectNormal;

    glm::vec3 interpolatedDiffuseColor;
    glm::vec3 interpolatedSpecularColor;
    float interpolatedShininess = 0.0f;
    float t = FLT_MAX;
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

inline bool isNearlyZero(float val)
{
    return abs(val) < eps;
}

#define SSAA_RANDOM 0
#define PRINT_INPUT 0
#endif //RAYTRACER_HEADER