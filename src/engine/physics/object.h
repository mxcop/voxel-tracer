#pragma once

/**
 * @brief Physics object data.
 */
struct PhyObject {
    /* Unique ID. */
    u32 uid;

    /* Physics attributes. */
    float3 position, velocity, force;
    f32 mass;

    PhyObject() = default;
    PhyObject(const float3 pos, const f32 mass) : position(pos), mass(mass) {}
};
