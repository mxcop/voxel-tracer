#pragma once

/* Collider types */
enum ColliderType { PLANE = 0, SPHERE = 1, BOX = 2 };

struct Collider {
    ColliderType type;

    Collider(ColliderType type) : type(type){};

    virtual float3 furthest_point(const Transform& t, const float3& dir) const { return t.position; };
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
