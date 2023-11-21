#include "Renderable.h"

// bool isNearlyZero(float val)
// {
//     return abs(val) < eps;
// }


class Triangle : public Renderable {
public:
    Vertex v[3];

    bool intersectRay(const glm::vec3& ray_o, const glm::vec3& ray_d, IntersectData& data) const 
    {
            // test if ray and plane are parallel
        glm::vec3 v1(this->v[0].position[0], this->v[0].position[1], this->v[0].position[2]);
        glm::vec3 v2(this->v[1].position[0], this->v[1].position[1], this->v[1].position[2]);
        glm::vec3 v3(this->v[2].position[0], this->v[2].position[1], this->v[2].position[2]);

        glm::vec3 planeNormal = glm::normalize(glm::cross(v2 - v1, v3 - v1));

        // if ray is parallel to the plane
        float nDotd = glm::dot(planeNormal, ray_d);
        if (abs(nDotd) < eps)
        {
            return false;
        }

        // calculate ray plane intersection
        float planeCoefficient_d = -(glm::dot(planeNormal, v1));
        float t = -(glm::dot(planeNormal, ray_o) + planeCoefficient_d) / nDotd;

        // intersection behind ray origin
        if (t <= eps)
        {
            return false;
        }

        // calculate plane intersection point
        glm::vec3 I = ray_o + t * ray_d;

        //try to 2D projection onto different planes and calculate barycentric coordinates
        float areaV1V2V3 = -1.0f;
        float areaV1V2I = -1.0f;
        float areaV2V3I = -1.0f;
        float areaV1IV3 = -1.0f;

        float weightV1V2I = -1.0f;
        float weightV2V3I = -1.0f;
        float weightV1IV3 = -1.0f;

        // try xy plane
        if (!isNearlyZero(glm::dot(planeNormal, glm::vec3(0.0f, 0.0f, 1.0f))))
        {
            areaV1V2V3 = 0.5f * ((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
            areaV1V2I = 0.5f * ((v2.x - v1.x) * (I.y - v1.y) - (I.x - v1.x) * (v2.y - v1.y));
            areaV2V3I = 0.5f * ((v2.x - I.x) * (v3.y - I.y) - (v3.x - I.x) * (v2.y - I.y));
            areaV1IV3 = 0.5f * ((I.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (I.y - v1.y));
        }
        // try xz plane
        else if (!isNearlyZero(glm::dot(planeNormal, glm::vec3(0.0f, 1.0f, 0.0f))))
        {
            areaV1V2V3 = 0.5f * ((v2.x - v1.x) * (v3.z - v1.z) - (v3.x - v1.x) * (v2.z - v1.z));
            areaV1V2I = 0.5f * ((v2.x - v1.x) * (I.z - v1.z) - (I.x - v1.x) * (v2.z - v1.z));
            areaV2V3I = 0.5f * ((v2.x - I.x) * (v3.z - I.z) - (v3.x - I.x) * (v2.z - I.z));
            areaV1IV3 = 0.5f * ((I.x - v1.x) * (v3.z - v1.z) - (v3.x - v1.x) * (I.z - v1.z));
        }
        // has to be yz plane
        else
        {
            areaV1V2V3 = 0.5f * ((v2.y - v1.y) * (v3.z - v1.z) - (v3.y - v1.y) * (v2.z - v1.z));
            areaV1V2I = 0.5f * ((v2.y - v1.y) * (I.z - v1.z) - (I.y - v1.y) * (v2.z - v1.z));
            areaV2V3I = 0.5f * ((v2.y - I.y) * (v3.z - I.z) - (v3.y - I.y) * (v2.z - I.z));
            areaV1IV3 = 0.5f * ((I.y - v1.y) * (v3.z - v1.z) - (v3.y - v1.y) * (I.z - v1.z));
        }

        // test if the intersection is in triangle, reject if we get negative weight
        weightV1V2I = areaV1V2I / areaV1V2V3;
        if (weightV1V2I < 0.0f)
        {
            return false;
        }

        weightV2V3I = areaV2V3I / areaV1V2V3;
        if (weightV2V3I < 0.0f)
        {
            return false;
        }

        weightV1IV3 = areaV1IV3 / areaV1V2V3;
        if (weightV1IV3 < 0.0f)
        {
            return false;
        }

        // fill out intersection data
        glm::vec3 v1_normal(this->v[0].normal[0], this->v[0].normal[1], this->v[0].normal[2]);
        glm::vec3 v2_normal(this->v[1].normal[0], this->v[1].normal[1], this->v[1].normal[2]);
        glm::vec3 v3_normal(this->v[2].normal[0], this->v[2].normal[1], this->v[2].normal[2]);

        glm::vec3 v1_diffuse(this->v[0].color_diffuse[0], this->v[0].color_diffuse[1], this->v[0].color_diffuse[2]);
        glm::vec3 v2_diffuse(this->v[1].color_diffuse[0], this->v[1].color_diffuse[1], this->v[1].color_diffuse[2]);
        glm::vec3 v3_diffuse(this->v[2].color_diffuse[0], this->v[2].color_diffuse[1], this->v[2].color_diffuse[2]);

        glm::vec3 v1_specular(this->v[0].color_specular[0], this->v[0].color_specular[1], this->v[0].color_specular[2]);
        glm::vec3 v2_specular(this->v[1].color_specular[0], this->v[1].color_specular[1], this->v[1].color_specular[2]);
        glm::vec3 v3_specular(this->v[2].color_specular[0], this->v[2].color_specular[1], this->v[2].color_specular[2]);

        float v1_shininess = this->v[0].shininess;
        float v2_shininess = this->v[1].shininess;
        float v3_shininess = this->v[2].shininess;

        // interpolate
        data.intersectPoint = I;
        data.intersectNormal = glm::normalize(v1_normal * weightV2V3I + v2_normal * weightV1IV3 + v3_normal * weightV1V2I);
        data.interpolatedDiffuseColor = v1_diffuse * weightV2V3I + v2_diffuse * weightV1IV3 + v3_diffuse * weightV1V2I;
        data.interpolatedSpecularColor = v1_specular * weightV2V3I + v2_specular * weightV1IV3 + v3_specular * weightV1V2I;
        data.interpolatedShininess = v1_shininess * weightV2V3I + v2_shininess * weightV1IV3 + v3_shininess * weightV1V2I;
        data.t = t;

        return true;
    
    }
};