#pragma once

class Pyramid {
    struct Plane {
        float3 normal;
        f32 d;
    };

    float3 origin = 0;
    Plane planes[4] = {};

   public:
    Pyramid() = default;
    /**
     * @brief Create a pyramid from an origin and 3 corner rays.
     *
     * @param o Frustum origin point.
     * @param tl Frustum direction top left.
     * @param tr Frustum direction top right.
     * @param bl Frustum direction bottom left.
     */
    Pyramid(const float3& o, float3 tl, float3 tr, float3 bl);

    float2 project(const float3& point) const;
};
