#pragma once

#include "graphics/primitives/traceable.h"

/**
 * @brief Capsule... not much to say here.
 */
struct Capsule : public Traceable {
    float3 a = 0, b = 0;
    f32 radius = 1;

    Capsule() = default;
    Capsule(const float3 a, const float3 b, const f32 radius);

    /* Trace-able functions */
    AABB get_aabb() const override;
    float3 center() const override;
    HitInfo intersect(const Ray& ray) const override;
};
