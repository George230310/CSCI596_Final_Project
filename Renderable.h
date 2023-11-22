#ifndef RENDERABLE_H
#define RENDERABLE_H
#include "raytracer.h"
class Renderable
{
    public:
        virtual bool intersectRay(const glm::vec3& ray_o, const glm::vec3& ray_d, IntersectData& data) const = 0;

        bool intersectSegment(const glm::vec3& a, const glm::vec3& b)
        {
            IntersectData data;
            glm::vec3 ab = b - a;
            glm::vec3 ray_d = glm::normalize(ab);

            if (this->intersectRay(a, ray_d, data))//(rayIntersectTriangle(a, ray_d, triangle, data))
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
};
#endif