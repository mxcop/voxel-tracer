#include "precomp.h"
#include "collision.h"

CollisionPoints::CollisionPoints(const float3& a, const float3& b) : a(a), b(b) {
    const float3 extend = b - a;
    dist = length(extend);
    normal = extend / dist;
    collision = true;
}

CollisionPoints::CollisionPoints(const float3& a, const float3& b, const float3& normal,
                                 const f32 dist)
    : a(a), b(b), normal(normal), dist(dist), collision(true) {}

CollisionPoints ct_sphere_vs_sphere(const Collider* a, const Transform* ta, const Collider* b,
                                    const Transform* tb) {
    using Sphere = SphereCollider;
    Sphere* A = (Sphere*)a;
    Sphere* B = (Sphere*)b;

    const float3 aCenter = A->center + ta->position;
    const float3 bCenter = B->center + tb->position;

    const float3 ab = bCenter - aCenter;

    const f32 a_radius = A->radius * major(ta->scale);
    const f32 b_radius = B->radius * major(tb->scale);

    const f32 distance = length(ab);

    if (distance < 0.00001f || distance > a_radius + b_radius) {
        return CollisionPoints();
    }

    const float3 normal = normalize(ab);

    const float3 a_deep = aCenter + normal * a_radius;
    const float3 b_deep = bCenter - normal * b_radius;

    return CollisionPoints(a_deep, b_deep);
}

CollisionPoints ct_plane_vs_sphere(const Collider* a, const Transform* ta, const Collider* b,
                                   const Transform* tb) {
    using Plane = PlaneCollider;
    using Sphere = SphereCollider;

    Plane* A = (Plane*)a;
    Sphere* B = (Sphere*)b;

    const float3 b_center = B->center + tb->position;
    const f32 b_radius = B->radius * major(tb->scale);

    const float3 normal = ta->rotation.rotateVector(A->normal);
    const float3 plane_d = normal * A->d + ta->position;

    /* Distance from center of sphere to plane surface */
    const f32 distance = dot(b_center - plane_d, normal);
    // const f32 d = dot(normal, b_center) - plane_d;

    if (distance > b_radius) {
        return CollisionPoints();
    }

    const float3 a_deep = b_center - normal * b_radius;
    const float3 b_deep = b_center - normal * distance;

    return CollisionPoints(a_deep, b_deep, normal, distance);
}
