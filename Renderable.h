#include "raytracer.h"
class Renderable
{
    public:
        virtual bool intersectRay(const glm::vec3& ray_o, const glm::vec3& ray_d, IntersectData& data) const = 0;

};