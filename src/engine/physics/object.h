#pragma once

#include "transform.h"
#include "collision/colliders.h"

/**
 * @brief Physics object data.
 */
struct PhyObject {
    /* Unique ID. */
    u32 uid = 0;

    /* Physics attributes. */
    float3 velocity = 0, force = 0;
    f32 mass = 0;

    Collider* collider = nullptr;
    Transform transform;

    PhyObject() = default;
    PhyObject(Collider* collider, const float3 pos, const f32 mass)
        : collider(collider), transform(Transform(pos)), mass(mass) {}
};
