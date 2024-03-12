#pragma once

/* Describes an objects location */
struct Transform {
    float3 position, scale;
    quat rotation;

    Transform() = default;
    explicit Transform(const float3& position) : position(position), scale(1) {}
};
