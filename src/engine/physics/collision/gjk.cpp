#include "gjk.h"

#include "simplex.h"

bool next_simplex(Simplex& points, float3& dir);

bool gjk(const Collider& a, const Transform& ta, const Collider& b, const Transform& tb) {
    float3 support = gjk_support(a, ta, b, tb, float3(1, 0, 0));

    /* Array of points, maximum count is 4 */
    Simplex points;
    points.push_front(support);

    /* Next direction towards the origin */
    float3 dir = -support;

    for (;;) {
        support = gjk_support(a, ta, b, tb, dir);

        if (dot(support, dir) <= 0) {
            return false; /* no collision */
        }

        points.push_front(support);

        if (next_simplex(points, dir)) {
            return true;
        }
    }
}

static inline bool same_dir(const float3& dir, const float3& ao) { return dot(dir, ao) > 0; }

/* Line case */
static bool simplex_line(Simplex& points, float3& dir) {
    const float3 a = points[0], b = points[1];
    const float3 ab = b - a, ao = -a;

    if (same_dir(ab, ao)) {
        dir = cross(cross(ab, ao), ab);
    } else {
        points = {a};
        dir = ao;
    }

    return false;
}

/* Triangle case */
static bool simplex_tri(Simplex& points, float3& dir) {
    const float3 a = points[0];
    const float3 b = points[1];
    const float3 c = points[2];

    const float3 ab = b - a;
    const float3 ac = c - a;
    const float3 ao = -a;

    const float3 abc = cross(ab, ac);

    if (same_dir(cross(abc, ac), ao)) {
        if (same_dir(ac, ao)) {
            points = {a, c};
            dir = cross(cross(ac, ao), ac);
        } else {
            return simplex_line(points = {a, b}, dir);
        }
    }

    else {
        if (same_dir(cross(ab, abc), ao)) {
            return simplex_line(points = {a, b}, dir);
        } else {
            if (same_dir(abc, ao)) {
                dir = abc;
            } else {
                points = {a, c, b};
                dir = -abc;
            }
        }
    }

    return false;
}

/* Tetrahedron case */
static bool simplex_tetra(Simplex& points, float3& dir) {
    const float3 a = points[0];
    const float3 b = points[1];
    const float3 c = points[2];
    const float3 d = points[3];

    const float3 ab = b - a;
    const float3 ac = c - a;
    const float3 ad = d - a;
    const float3 ao = -a;

    const float3 abc = cross(ab, ac);
    const float3 acd = cross(ac, ad);
    const float3 adb = cross(ad, ab);

    if (same_dir(abc, ao)) {
        return simplex_tri(points = {a, b, c}, dir);
    }

    if (same_dir(acd, ao)) {
        return simplex_tri(points = {a, c, d}, dir);
    }

    if (same_dir(adb, ao)) {
        return simplex_tri(points = {a, d, b}, dir);
    }

    return true;
}

bool next_simplex(Simplex& points, float3& dir) {
    switch (points.get_size()) {
        case 2:
            return simplex_line(points, dir);
        case 3:
            return simplex_tri(points, dir);
        case 4:
            return simplex_tetra(points, dir);
    }

    return false; /* ! unreachable ! */
}
