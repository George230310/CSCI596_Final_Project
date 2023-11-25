#ifndef RAYTRACER_HEADER
#define RAYTRACER_HEADER

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include "glm/glm.hpp"

#include "png++/png.hpp"
#include <math.h>
#include <vector>
#include <random>

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

#define SSAAA_RANDOM 1
#endif //RAYTRACER_HEADER