#pragma once

#include "colliders.h"
#include "../transform.h"

/**
 * @brief Represents a potential collision between two primitives as two points.
 */
struct CollisionPoints {
    float3 a = 0;      /* Furthest point of A into B */
    float3 b = 0;      /* Furthest point of B into A */
    float3 normal = 0; /* Direction from A to B */
    float dist = 0;    /* Distance between A and B */
    bool collision = false;

    CollisionPoints() = default;
    CollisionPoints(const float3& a, const float3& b);
    CollisionPoints(const float3& a, const float3& b, const float3& normal, const f32 dist);
};

/**
 * @brief Represents a collision between two objects.
 */
struct Collision {
    PhyObject* a;
    PhyObject* b;
    CollisionPoints points;

    Collision(PhyObject* a, PhyObject* b, const CollisionPoints points)
        : a(a), b(b), points(points){};
};

/*
 * ===== Collision Tests =====
 *      "ct_{p1}_vs_{p2}"
 */

CollisionPoints ct_sphere_vs_sphere(const Collider* a, const Transform* ta, const Collider* b,
                                    const Transform* tb);

CollisionPoints ct_plane_vs_sphere(const Collider* a, const Transform* ta, const Collider* b,
                                   const Transform* tb);

CollisionPoints ct_box_vs_sphere(const Collider* a, const Transform* ta, const Collider* b,
                                 const Transform* tb);

/* Function ptr */
using FindContactFunc = CollisionPoints (*)(const Collider*, const Transform*, const Collider*,
                                            const Transform*);

const FindContactFunc func_table[3][3] = {
    /* Plane, Sphere, Box */
    {nullptr, nullptr, nullptr},                                 /* Plane */
    {ct_sphere_vs_sphere, ct_plane_vs_sphere, ct_box_vs_sphere}, /* Sphere */
    {nullptr, nullptr, nullptr}                                  /* Box */
};

inline CollisionPoints collision_test(const Collider* a, const Transform* at, const Collider* b,
                                      const Transform* bt) {
    /* Swap "Plane vs Sphere" into "Sphere vs Plane" */
    bool swap = b->type < a->type;
    if (swap) {
        std::swap(a, b);
        std::swap(at, bt);
    }

    /* Dispatch the collision test */
    CollisionPoints points = func_table[a->type][b->type](a, at, b, bt);

    /* Swap the results back if the input was swapped */
    if (swap) {
        std::swap(points.a, points.b);
        points.normal = -points.normal;
    }

    return points;
}
