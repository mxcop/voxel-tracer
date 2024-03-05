#pragma once

#include "graphics/primitives/traceable.h"

/**
 * @brief Sphere... not much to say here.
 */
struct Sphere : public Traceable {
    /* Center point */
    float3 pos = 0;
    f32 radius = 1;

    Sphere() = default;
    Sphere(const float3 pos, const f32 radius);

    /* Trace-able functions */
    AABB get_aabb() const override;
    float3 center() const override;
    HitInfo intersect(const Ray& ray) const override;
};
