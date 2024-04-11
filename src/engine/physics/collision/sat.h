#pragma once

/**
 * @brief Project a 3D point onto a position and unit vector.
 *
 * @param p The point to project.
 * @param o The origin of the unit vector.
 * @param n The unit vector. (MUST BE A UNIT VECTOR)
 * @return The projection of P onto N.
 */
inline f32 project_onto(const float3& p, const float3& o, const float3& n) {
    const float3 local_p = p - o;
    return dot(n, local_p);
}

/**
 * @brief Project a 3D point onto a vector. (good enough for comparisons)
 *
 * @param p The point to project.
 * @param n The unit vector to project onto.
 * @return The projection of P onto N.
 */
inline f32 project_onto(const float3& p, const float3& n) { return dot(p, n); }

struct box_t {
    float3 min, max;
    quat rot;
};

struct corners {
    float3& operator[](const i32 i) { return corners[i]; }
    const float3& operator[](const i32 i) const { return corners[i]; }

   private:
    float3 corners[8];
};

/**
 * @brief Get the 8 corners of a box.
 */
inline corners box_corners(const box_t& box) {
    corners r;
    for (u32 c = 0; c < 8; ++c) {
        /* Get a corner */
        r[c].x = (c & 0b001) ? box.max.x : box.min.x;
        r[c].y = (c & 0b010) ? box.max.y : box.min.y;
        r[c].z = (c & 0b100) ? box.max.z : box.min.z;
        /* WARN: no pivot used here... */
        // r[c] = box.rot.rotate_vec(r[c]);
    }
    return r;
}

/**
 * @brief Find the min and max value of the projected box.
 */
inline void projected_minmax(const corners& box, const float3& n, f32& min, f32& max) {
    min = BIG_F32, max = -BIG_F32;
    for (u32 c = 0; c < 8; ++c) {
        /* Save min and max projection */
        const f32 proj = project_onto(box[c], n);
        min = fminf(min, proj), max = fmaxf(max, proj);
    }
}

/**
 * @brief Find if there is a separating plane between two sets of corners.
 */
inline bool separation(const corners& corners_a, const corners& corners_b, const float3& n) {
    /* Get min and max projection of A and B onto the normal */
    f32 a_min, a_max, b_min, b_max;
    projected_minmax(corners_a, n, a_min, a_max);
    projected_minmax(corners_b, n, b_min, b_max);

    /* Check for separation */
    if (a_max >= b_min && b_max >= a_min) {
        return false;
    }
    return true; /* Separating plane found! */
}

/**
 * @brief Find if two boxes are intersecting each other.
 */
inline bool box_sat(const box_t& box_a, const box_t& box_b) {
    /* Look for a separation on all 15 axes */
    const corners corners_a = box_corners(box_a);
    const corners corners_b = box_corners(box_b);

    /* Check box A normals (3 axes) */
    const quat_axes rot_a = box_a.rot.axes();
    for (u32 n = 0; n < 3; ++n) {
        /* Get a rotated normal */
        const float3 normal = rot_a[n];

        /* Check for separation */
        if (separation(corners_a, corners_b, normal)) {
            return false; /* Separating plane found! */
        }
    }

    /* Check box B normals (3 axes) */
    const quat_axes rot_b = box_b.rot.axes();
    for (u32 n = 0; n < 3; ++n) {
        /* Get a rotated normal */
        const float3 normal = rot_b[n];

        /* Check for separation */
        if (separation(corners_a, corners_b, normal)) {
            return false; /* Separating plane found! */
        }
    }

    /* Check cross products (9 axes) */
    for (u32 an = 0; an < 3; ++an) {
        /* Get a rotated normal */
        const float3 normal_a = rot_a[an];

        for (u32 bn = 0; bn < 3; ++bn) {
            /* Get a rotated normal */
            const float3 normal_b = rot_b[bn];

            /* Cross projection normal */
            const float3 normal = normalize(cross(normal_a, normal_b));

            /* Check for separation */
            if (separation(corners_a, corners_b, normal)) {
                return false; /* Separating plane found! */
            }
        }
    }

    return true;
}

/**
 * @brief Find if there is a separating plane between two sets of corners.
 */
inline bool separation(const corners& corners_a, const Pyramid& pyramid, const float3& n) {
    /* Get min and max projection of A and B onto the normal */
    f32 a_min, a_max, b_min, b_max;
    projected_minmax(corners_a, n, a_min, a_max);
    pyramid.projected_minmax(n, b_min, b_max);

    /* Check for separation */
    if (a_max >= b_min && b_max >= a_min) {
        return false;
    }
    return true; /* Separating plane found! */
}

inline f32 projected_dist(const corners& corners_a, const Pyramid& pyramid) {
    const float3 n = pyramid.forward;
    f32 a_min, a_max;
    projected_minmax(corners_a, n, a_min, a_max);
    const f32 proj = project_onto(pyramid.origin, n);
    const f32 dist = a_min - proj;
    return dist;
}

/**
 * @brief Find if a box and a pyramid are intersecting each other.
 */
inline f32 box_pyramid_sat(const box_t& box, const Pyramid& pyramid) {
    /* Look for a separation on all 25 axes */
    const corners corners_a = box_corners(box);

    /* Check pyramid normals (4 axes) */
    for (u32 n = 0; n < 4; ++n) {
        /* Get a rotated normal */
        const float3 normal = pyramid.planes[n].normal;

        /* Check for separation */
        if (separation(corners_a, pyramid, normal)) {
            return BIG_F32; /* Separating plane found! */
        }
    }

    /* Check box normals (3 axes) */
    const quat_axes rot_a = box.rot.axes();
    for (u32 n = 0; n < 3; ++n) {
        /* Get a rotated normal */
        const float3 normal = rot_a[n];

        /* Check for separation */
        if (separation(corners_a, pyramid, normal)) {
            return BIG_F32; /* Separating plane found! */
        }
    }

/* Check cross products (18 axes) */
#if ACCURATE_PYRAMID_TRACING
    for (u32 an = 0; an < 3; ++an) {
        /* Get a rotated normal */
        const float3 normal_a = rot_a[an];

        for (u32 bn = 0; bn < 6; ++bn) {
            /* Get a rotated normal */
            const float3 normal_b = pyramid.rays[bn];

            /* Cross projection normal */
            const float3 normal = normalize(cross(normal_b, normal_a));

            /* Check for separation */
            if (separation(corners_a, pyramid, normal)) {
                return BIG_F32; /* Separating plane found! */
            }
        }
    }
#endif

    return projected_dist(corners_a, pyramid);
}
