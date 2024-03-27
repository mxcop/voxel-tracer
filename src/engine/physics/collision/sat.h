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
inline f32 project_onto(const float3& p, const float3& n) { return dot(n, p); }

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
        r[c] = box.rot.rotate_vec(r[c]);
    }
    return r;
}

/**
 * @brief Find the min and max value of the projected box.
 */
inline void projected_minmax(const corners& box, const float3& n, f32& min, f32& max) {
    for (u32 c = 0; c < 8; ++c) {
        /* Save min and max projection */
        const f32 proj = project_onto(box[c], n);
        min = fminf(min, proj), max = fmaxf(max, proj);
    }
}

/**
 * @brief Find of two boxes are intersecting each other.
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

        /* Get min and max projection of A and B onto the normal */
        f32 a_min, a_max, b_min, b_max;
        projected_minmax(corners_a, normal, a_min, a_max);
        projected_minmax(corners_b, normal, b_min, b_max);

        /* Check for separation */
        if (a_max < b_min || b_max < a_min) {
            return false; /* Separating plane found! */
        }
    }

    /* Check box B normals (3 axes) */
    const quat_axes rot_b = box_b.rot.axes();
    for (u32 n = 0; n < 3; ++n) {
        /* Get a rotated normal */
        const float3 normal = rot_b[n];

        /* Get min and max projection of A and B onto the normal */
        f32 a_min, a_max, b_min, b_max;
        projected_minmax(corners_a, normal, a_min, a_max);
        projected_minmax(corners_b, normal, b_min, b_max);

        /* Check for separation */
        if (a_max < b_min || b_max < a_min) {
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

            /* Get min and max projection of A and B onto the normal */
            f32 a_min, a_max, b_min, b_max;
            projected_minmax(corners_a, normal, a_min, a_max);
            projected_minmax(corners_b, normal, b_min, b_max);

            /* Check for separation */
            if (a_max < b_min || b_max < a_min) {
                return false; /* Separating plane found! */
            }
        }
    }

    return true;
}
