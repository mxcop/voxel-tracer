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
bool box_sat(const box_t& box_a, const box_t& box_b) {
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
            const float3 normal = cross(normal_a, normal_b);

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
inline bool separation(const corners& corners_a, const Pyramid& pyramid, const float3& n, bool debug = false) {
    /* Get min and max projection of A and B onto the normal */
    f32 a_min, a_max, b_min, b_max;
    projected_minmax(corners_a, n, a_min, a_max);
    pyramid.projected_minmax(n, b_min, b_max); // TODO: this projection seems wrong...

    if (debug) {
        const float3 center = corners_a[0] + (corners_a[7] - corners_a[0]) * 0.5f;
        //const float3 center = {5,0,5};
        db::draw_line(a_min * n, a_max * n, 0xFFFF0000);
        db::draw_line(b_min * n, b_max * n, 0xFF0000FF);
    }

    /* Check for separation */
    if (a_max >= b_min && b_max >= a_min) {
        return false;
    }
    return true; /* Separating plane found! */
}

/**
 * @brief Find if a box and a pyramid are intersecting each other.
 */
bool box_pyramid_sat(const box_t& box, const Pyramid& pyramid, bool debug = false) {
    /* Look for a separation on all 15 axes */
    const corners corners_a = box_corners(box);
    
    /* Check box normals (3 axes) */
    const quat_axes rot_a = box.rot.axes();
    for (u32 n = 0; n < 3; ++n) {
        /* Get a rotated normal */
        const float3 normal = rot_a[n];

        if (debug) db::draw_normal(box.min + (box.max - box.min) * 0.5f, normal, 0xFFFF0000);

        /* Check for separation */
        if (separation(corners_a, pyramid, normal)) {
            return false; /* Separating plane found! */
        }
    }

    /* Check pyramid normals (5 axes) */
    if (separation(corners_a, pyramid, pyramid.forward)) {
        return false; /* Separating plane found! */
    }

    for (u32 n = 0; n < 4; ++n) {
        /* Get a rotated normal */
        const float3 normal = pyramid.planes[n].normal;

        if (debug) db::draw_normal(box.min + (box.max - box.min) * 0.5f, normal, 0xFF0000FF);

        /* Check for separation */
        if (separation(corners_a, pyramid, normal)) {
            return false; /* Separating plane found! */
        }
    }

    //if (separation(corners_a, pyramid, (pyramid.planes[0].normal + pyramid.planes[2].normal) * 0.5f,
    //               debug)) {
    //    return false; /* Separating plane found! */
    //}
    //if (separation(corners_a, pyramid, (pyramid.planes[1].normal + pyramid.planes[2].normal) * 0.5f,
    //               debug)) {
    //    return false; /* Separating plane found! */
    //}
    //if (separation(corners_a, pyramid, (pyramid.planes[0].normal + pyramid.planes[3].normal) * 0.5f,
    //               debug)) {
    //    return false; /* Separating plane found! */
    //}
    //if (separation(corners_a, pyramid, (pyramid.planes[1].normal + pyramid.planes[3].normal) * 0.5f,
    //               debug)) {
    //    return false; /* Separating plane found! */
    //}

    /* Check cross products (12 axes) */
    for (u32 an = 0; an < 3; ++an) {
        /* Get a rotated normal */
        const float3 normal_a = rot_a[an];

        /* Forward normal axis */
        //const float3 f_normal = normalize(cross(pyramid.forward, normal_a));

        //if (debug) db::draw_normal(box.min + (box.max - box.min) * 0.5f, f_normal, 0xFFFF00FF);

        //if (separation(corners_a, pyramid, f_normal, debug)) {
        //    return false; /* Separating plane found! */
        //}

        for (u32 bn = 0; bn < 6; ++bn) {
            /* Get a rotated normal */
            const float3 normal_b = pyramid.rays[bn];

            /* Cross projection normal */
            const float3 normal = normalize(cross(normal_b, normal_a));

            if (debug) db::draw_normal(box.min + (box.max - box.min) * 0.5f, normal, 0xFFFF00FF);

            /* Check for separation */
            if (separation(corners_a, pyramid, normal)) {
                return false; /* Separating plane found! */
            }
        }
    }

    return true;
}
