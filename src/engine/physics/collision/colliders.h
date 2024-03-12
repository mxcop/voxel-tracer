#pragma once

/* Collider types */
enum ColliderType { PLANE = 0, SPHERE = 1 };

struct Collider {
    ColliderType type;

    Collider(ColliderType type) : type(type){};
};

/*
 * ===== Collider Definitions =====
 */

struct SphereCollider : Collider {
    float3 center;
    f32 radius;

    SphereCollider() = delete;
    SphereCollider(const float3& center, const f32 radius) : Collider(SPHERE), center(center), radius(radius) {};
};

struct PlaneCollider : Collider {
    float3 normal;
    f32 d;

    PlaneCollider() = delete;
    PlaneCollider(const float3& normal, const f32 d) : Collider(PLANE), normal(normal), d(d){};
};
