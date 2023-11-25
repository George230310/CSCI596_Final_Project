#include "Renderable.h"

class Sphere : public Renderable {
public:
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
    Sphere ()
    {
        
    }
    Sphere(double position[3], double color_diffuse[3], double color_specular[3], double shininess, double radius);
    

    // Implementation of the intersectRay function
    bool intersectRay(const glm::vec3& ray_o, const glm::vec3& ray_d, IntersectData& data) const
    {
        glm::vec3 sphereCenter(this->position[0], this->position[1], this->position[2]);
        // if the radius of sphere is zero, no intersection
        if (this->radius < eps)
        {
            return false;
        }

        // float b = 2.0f * (ray_d.x * (ray_o.x - this->position[0]) + ray_d.y * (ray_o.y - this->position[1]) + ray_d.z * (ray_o.z - this->position[2]));
        float b = 2.0f * dot(ray_d,(ray_o - sphereCenter));

        // float c = powf(ray_o.x - this->position[0], 2.0f) + powf(ray_o.y - this->position[1], 2.0f) + powf(ray_o.z - this->position[2], 2.0f) - powf(this->radius, 2.0f);
        float c = powf(length(ray_o-sphereCenter),2.0f) - powf(this->radius, 2.0f);

        float delta = powf(b, 2.0f) - 4.0f * c;

        if (delta < 0.0f)
        {
            return false;
        }

        float t0 = (-b + sqrtf(delta)) / 2.0f;
        float t1 = (-b - sqrtf(delta)) / 2.0f;

        // if both t are less or equal to zero, no intersection
        if (t0 < eps && t1 < eps)
        {
            return false;
        }

        float t_intersect = 0.0f;

        if (t0 < eps)
        {
            t_intersect = t1;
        }
        else if (t1 < eps)
        {
            t_intersect = t0;
        }
        else
        {
            t_intersect = std::min(t0, t1);
        }

        // fill out intersection data
        data.intersectPoint = ray_o + ray_d * t_intersect;
        
        data.intersectNormal = (data.intersectPoint - sphereCenter) / (float)this->radius;

        // negate if ray originates inside the sphere
        if (glm::length(ray_o - sphereCenter) < (float)this->radius)
        {
            data.intersectNormal = -data.intersectNormal;
        }

        data.t = t_intersect;

        data.interpolatedDiffuseColor = glm::vec3(this->color_diffuse[0], this->color_diffuse[1], this->color_diffuse[2]);
        data.interpolatedSpecularColor = glm::vec3(this->color_specular[0], this->color_specular[1], this->color_specular[2]);
        data.interpolatedShininess = (float)this->shininess;

        return true;
    }
};
