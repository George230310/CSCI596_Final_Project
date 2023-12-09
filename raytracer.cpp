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
#include "libs/glm/gtc/random.hpp"
#include "mpi.h"

char * filename = NULL;
int mode = MODE_DISPLAY;
unsigned int buffer[HEIGHT][WIDTH][3];

// parallelization
int myId;
int nProc;

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

// function to render a rectangular block of the super scaled image from indicated start column & row to end column & row
void renderBlock(int startCol, int startRow, int endCol, int endRow)
{
    int rayIndex = (startCol + 1) * (startRow + 1) - 1;

    // go through each pixel
    for (int i = startCol; i < endCol; ++i)
    {
        for (int j = startRow; j < endRow; ++j)
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
}

// Parallelization: partition the scene into horizontal strips

void processScene_manager()
{
    double renderTime = 0.0f;
    double startTime;
    double finishTime;

    startTime = MPI_Wtime();

    // sub-divide light sources
    SubdivideLightSources();

    // process the scene
    MPI_Status status;
    int scaledImageHeight = HEIGHT * SSAA_Coefficient;
    int scaledImageWidth = WIDTH * SSAA_Coefficient;

    // buffer to receive pixels from other processes
    int bufferItemCount = 3 * scaledImageHeight * scaledImageWidth;
    int* recvPixels = new int[bufferItemCount];

    // TODO: the manager renders the first strip and any left overs, probably better way to do this
    int rowsPerBlock = scaledImageHeight / nProc;
    int leftOverRows = scaledImageHeight % nProc;
    renderBlock(0, 0, scaledImageWidth, rowsPerBlock);

    if(leftOverRows > 0)
    {
        renderBlock(0, scaledImageHeight - leftOverRows + 2, scaledImageWidth, scaledImageHeight);
    }

    // receive from workers
    for(int i = 1; i < nProc; i++)
    {
        MPI_Recv(recvPixels, bufferItemCount, MPI_INT, i, i, MPI_COMM_WORLD, &status);
        int rowOffset = i * rowsPerBlock;
        for(int a = 0; a < rowsPerBlock; a++)
        {
            for(int b = 0; b < scaledImageWidth; b++)
            {
                int row = a + rowOffset;
                int index = 3 * (row * scaledImageWidth + b);
                superScaledAllPixels[b][a][0] = recvPixels[index];
                superScaledAllPixels[b][a][1] = recvPixels[index + 1];
                superScaledAllPixels[b][a][2] = recvPixels[index + 2];
            }
        }
    }

    // TODO: manager solely responsible for scale down now, but there's definitely better way
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

    finishTime = MPI_Wtime();
    renderTime = finishTime - startTime;
    std::cout << "Execution Time: " << renderTime << " seconds" << std::endl;
}

void processScene_worker()
{
    // buffer to send pixels from other processes
    int scaledImageHeight = HEIGHT * SSAA_Coefficient;
    int scaledImageWidth = WIDTH * SSAA_Coefficient;
    int bufferItemCount = 3 * scaledImageHeight * scaledImageWidth;
    int* sendPixels = new int[bufferItemCount];

    // render my part
    int rowsPerBlock = scaledImageHeight / nProc;
    int rowOffset = myId * rowsPerBlock;
    renderBlock(0, rowOffset, scaledImageWidth, rowOffset + rowsPerBlock);

    // send the result to manager
    for(int a = 0; a < rowsPerBlock; a++)
    {
        for(int b = 0; b < scaledImageWidth; b++)
        {
            int row = a + rowOffset;
            int index = 3 * (row * scaledImageWidth + b);
            sendPixels[index] = superScaledAllPixels[b][a][0];
            sendPixels[index + 1] = superScaledAllPixels[b][a][1];
            sendPixels[index + 2] = superScaledAllPixels[b][a][2];
        }
    }

    MPI_Send(sendPixels, bufferItemCount, MPI_INT, 0, myId, MPI_COMM_WORLD);
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

  // output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      plot_pixel(x, y, allPixels[x][y][0], allPixels[x][y][1], allPixels[x][y][2]);
    }
  }
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

void save_png()
{
    std::vector<unsigned char> image;
    for(unsigned y = 0; y < HEIGHT; y++)
    {
        for(unsigned x = 0; x < WIDTH; x++)
        {
            image.push_back((unsigned char)buffer[y][x][0]);
            image.push_back((unsigned char)buffer[y][x][1]);
            image.push_back((unsigned char)buffer[y][x][2]);
            image.push_back(255);
        }
    }

    encodeOneStep(filename, image, WIDTH, HEIGHT);
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

int main(int argc, char ** argv)
{
  // MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  if (argc != 3)
  {  
    printf ("Usage: %s <input scenefile> <output pngname>\n", argv[0]);
    exit(0);
  }

  filename = argv[2];
  loadScene(argv[1]);

  // generate all rays from COP (should be run by all processes)
  generateAllRaysFromCOP();

  if(myId == 0)
  {
    // the manager process
    processScene_manager();

    // plot pixel to the buffer
    for(unsigned int x=0; x<WIDTH; x++)
    {
        for(unsigned int y=0; y<HEIGHT; y++)
        {
            plot_pixel(x, y, allPixels[x][y][0], allPixels[x][y][1], allPixels[x][y][2]);
        }
    }

    // save image from the buffer
    save_png();
  }
  else
  {
    // the worker processes
    // processScene_worker();
  }


  std::cout << "Proc #" << myId << " is done." << std::endl;
  MPI_Finalize();
  return 0;
}

