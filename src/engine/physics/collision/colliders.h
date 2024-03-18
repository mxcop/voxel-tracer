#pragma once

/* Collider types */
enum ColliderType { PLANE = 0, SPHERE = 1, BOX = 2, VOXELS = 3 };

struct Collider {
    ColliderType type;

    Collider(ColliderType type) : type(type){};

    virtual float3 furthest_point(const Transform& t, const float3& dir) const {
        return t.position;
    };
};

/*
 * ===== Collider Definitions =====
 */

struct SphereCollider : Collider {
    float3 center;
    f32 radius;

    SphereCollider() = delete;
    SphereCollider(const float3& center, const f32 radius)
        : Collider(SPHERE), center(center), radius(radius){};

    float3 furthest_point(const Transform& t, const float3& dir) const override;
};

struct PlaneCollider : Collider {
    float3 normal;
    f32 d;

    PlaneCollider() = delete;
    PlaneCollider(const float3& normal, const f32 d) : Collider(PLANE), normal(normal), d(d){};
};

struct BoxCollider : Collider {
    float3 min, max;

    BoxCollider() = delete;
    BoxCollider(const float3& min, const float3& max) : Collider(BOX), min(min), max(max){};

    float3 furthest_point(const Transform& t, const float3& dir) const override;
};

struct VoxelCollider : Collider {
    /* Size in world units */
    float3 extend;
    /* Size in voxels */
    int3 size;

    /* Abstracted voxel data */
    u8* voxels = nullptr;
    /* -> 2 bits per voxel? <-
     * 0 = Empty or surrounded. (never generates a contact point, and never needs to be checked)
     * 1 = Face voxel.          (never generates a contact point)
     * 2 = Edge voxel.          (can generate contact point, with "edge" or "corner")
     * 3 = Corner voxel.        (can generate contact point, with "face", "edge", or "corner")
     */

    VoxelCollider() = delete;
    VoxelCollider(const int3& grid_size, const f32 vpu, const u8* voxel_data);
};
