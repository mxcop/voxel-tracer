#pragma once

#include "graphics/primitives/traceable.h"

/**
 * @brief Oriented bounding box. (Box with rotation)
 */
struct OBB : public Traceable {
    /* Model matrix */
    mat4 model;
    float3 pos = 0, size = 0;

    OBB() = default;
    OBB(const float3 pos, const float3 size, const float3 axis = 0, const f32 angle = 0);

    /* Trace-able functions */
    AABB get_aabb() const override;
    float3 center() const override;
    HitInfo intersect(const Ray& ray) const override;

    f32 area() const;
};
