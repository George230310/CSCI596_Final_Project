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
#include "glm/gtc/random.hpp"

char * filename = NULL;
int mode = MODE_DISPLAY;
unsigned int buffer[HEIGHT][WIDTH][3];




struct Light
{
  double position[3];
  double color[3];
};

Renderable* renderables[MAX_SPHERES + MAX_TRIANGLES];
Light lights[MAX_LIGHTS];
std::vector<Light> subdividedLights;
double ambient_light[3];

// array of screen rays normalized directions
std::vector<glm::vec3> normalizedScreenRaysDirections;

// record pixel colors calculated
unsigned int allPixels[WIDTH][HEIGHT][3];

int num_triangles = 0;
int num_spheres = 0;
int num_renderables = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);


// ****************************** PROGRAM PARAMETERS *************************************
const unsigned int SSAA_Coefficient = 10;

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
        // go through all renderables in scene and find the intersection with the smallest t
        IntersectData data;
        IntersectData tempData;
        for (int k = 0; k < num_renderables; ++k)
        {
            if (renderables[k]->intersectRay(ray_o, ray_d, tempData)) 
            {
                // if we see a smaller t update intersection data
                if (tempData.t < data.t)
                {
                    data = tempData;
                }
            }
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

                // test whether we're shadowed by a renderable
                for (int n = 0; n < num_renderables; ++n)
                {
                    if (renderables[n]->intersectSegment(a,b))//(segmentIntersectSphere(a, b, renderables[n]))
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
            // go through all renderables in scene and find the intersection with the smallest t
            IntersectData data;
            IntersectData tempData;
            for (int k = 0; k < num_renderables; ++k)
            {
                if(renderables[k]->intersectRay(glm::vec3(0.0f, 0.0f, 0.0f), normalizedScreenRaysDirections[rayIndex], tempData))
                {
                    // if we see a smaller t update intersection data
                    if (tempData.t < data.t)
                    {
                        data = tempData;
                    }
                }
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

                    // test whether we're shadowed by a renderable
                    for (int n = 0; n < num_renderables; ++n)
                    {
                        if (renderables[n]->intersectSegment(a, b))
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

// function to generate Base ray from camera origin to pixel on virtual screen
glm::vec3 generateBaseRay(int i, int j)
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
    glm::vec3 pixelWidthOffset = (topRightScreenCorner - topLeftScreenCorner) / ((float)WIDTH);
    glm::vec3 pixelHeightOffset = (bottomRightScreenCorner - topRightScreenCorner) / ((float)HEIGHT);

    glm::vec3 ray = topLeftScreenCorner + (float)i * pixelWidthOffset + (float)j * pixelHeightOffset + pixelWidthOffset / 2.0f + pixelHeightOffset / 2.0f;
    ray = glm::normalize(ray);
    return ray;
}

// Anti-aliasing using random rays
void processScene_randomized()
{
    SubdivideLightSources();

    // Constants to project rays to each pixel on virtual screen
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
    glm::vec3 pixelWidthOffset = (topRightScreenCorner - topLeftScreenCorner) / ((float)WIDTH);
    glm::vec3 pixelHeightOffset = (bottomRightScreenCorner - topRightScreenCorner) / ((float)HEIGHT);

    // go through each pixel
    for (int i = 0; i < WIDTH ; ++i)
    {
        for (int j = 0; j < HEIGHT; ++j)
        {
            glm::vec3 baseRay = glm::normalize(topLeftScreenCorner + (float)i * pixelWidthOffset + (float)j * pixelHeightOffset + pixelWidthOffset / 2.0f + pixelHeightOffset / 2.0f);
            for(int l=0;l<SSAA_Coefficient*SSAA_Coefficient;l++)
            {
              glm::vec3 randomOffset  = glm::vec3(glm::gaussRand<float>(0.0f,0.008f), glm::gaussRand<float>(0.0f,0.008f), glm::gaussRand<float>(0.0f,0.008f));
              glm::vec3 finalRay = baseRay + randomOffset;
              IntersectData data;
              IntersectData tempData;
              // go through all renderables in scene and find the intersection with the smallest t
              for (int k = 0; k < num_renderables; ++k)
              {
                  if(renderables[k]->intersectRay(glm::vec3(0.0f, 0.0f, 0.0f), finalRay, tempData))
                  {
                      // if we see a smaller t update intersection data
                      if (tempData.t < data.t)
                      {
                          data = tempData;
                      }
                  }
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

                      // test whether we're shadowed by a renderable
                      for (int n = 0; n < num_renderables; ++n)
                      {
                          if (renderables[n]->intersectSegment(a, b))
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
                      allPixels[i][j][0] += (unsigned int)(255.0f * pixelColor.x);
                      allPixels[i][j][1] += (unsigned int)(255.0f * pixelColor.y);
                      allPixels[i][j][2] += (unsigned int)(255.0f * pixelColor.z);
                  }
                  else
                  {
                      glm::vec3 newDirection = glm::normalize(glm::reflect(finalRay, data.intersectNormal));
                      glm::vec3 reflectedColor = recursiveRayTrace(data.intersectPoint, newDirection, recursiveDepth);
                      glm::vec3 finalColor = (1.0f - recursiveReflection_lambda) * pixelColor + recursiveReflection_lambda * reflectedColor;

                      allPixels[i][j][0] += (unsigned int)(255.0f * finalColor.x);
                      allPixels[i][j][1] += (unsigned int)(255.0f * finalColor.y);
                      allPixels[i][j][2] += (unsigned int)(255.0f * finalColor.z);
                  }
                  
              }
              // else assign the background color
              else
              {
                  allPixels[i][j][0] += 255;
                  allPixels[i][j][1] += 255;
                  allPixels[i][j][2] += 255;
              }
            }
            allPixels[i][j][0]/= SSAA_Coefficient*SSAA_Coefficient;
            allPixels[i][j][1]/= SSAA_Coefficient*SSAA_Coefficient;
            allPixels[i][j][2]/= SSAA_Coefficient*SSAA_Coefficient;

        }
    }

  
}

//MODIFY THIS FUNCTION
void draw_scene()
{
    // generate all rays from COP
    generateAllRaysFromCOP();

    // process the scene
    if(SSAAA_RANDOM)
        processScene_randomized();
    else
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

void save_png()
{
    png::image< png::rgb_pixel > image(WIDTH, HEIGHT);
    for (png::uint_32 y = 0; y < image.get_height(); ++y)
    {
        for (png::uint_32 x = 0; x < image.get_width(); ++x)
        {
            image[HEIGHT-y-1][x] = png::rgb_pixel(buffer[y][x][0], buffer[y][x][1], buffer[y][x][2]);
            // non-checking equivalent of image.set_pixel(x, y, ...);
        }
    }
    image.write(filename);
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
void read_doubles(FILE* file, const char *check, double p[3])
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
      Triangle* newTriangle = new Triangle;
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",newTriangle->v[j].position);
        parse_doubles(file,"nor:",newTriangle->v[j].normal);
        parse_doubles(file,"dif:",newTriangle->v[j].color_diffuse);
        parse_doubles(file,"spe:",newTriangle->v[j].color_specular);
        parse_shi(file,&newTriangle->v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
        num_triangles++;
        renderables[num_renderables++] = newTriangle;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");
      Sphere* newSphere= new Sphere();
      
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
        num_spheres++;
        renderables[num_renderables++] = newSphere;
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

