#pragma once

/* SIMD (SSE) Ray packet structure. */
struct RayPacket {
    union {
        f128 ro[3]; /* Ray origin (XYZ) */
        struct {
            f128 ro_x;
            f128 ro_y;
            f128 ro_z;
        };
    };
    union {
        f128 rd[3]; /* Ray direction (XYZ) */
        struct {
            f128 rd_x;
            f128 rd_y;
            f128 rd_z;
        };
    };
    union {
        f128 ird[3]; /* Ray direction reciprocal (XYZ) */
        struct {
            f128 ird_x;
            f128 ird_y;
            f128 ird_z;
        };
    };
    float3 signs; /* Ray signs are shared by all rays in packet! */

    RayPacket() = default;
    inline RayPacket(const f128 shared_origin[3], const f128 dir[3]) {
        ro[0] = shared_origin[0];
        ro[1] = shared_origin[1];
        ro[2] = shared_origin[2];
        rd[0] = dir[0];
        rd[1] = dir[1];
        rd[2] = dir[2];
        ird[0] = _mm_rcp_ps(dir[0]);
        ird[1] = _mm_rcp_ps(dir[1]);
        ird[2] = _mm_rcp_ps(dir[2]);
        signs.x = copysign(1.0f, _mm_cvtss_f32(dir[0]));
        signs.y = copysign(1.0f, _mm_cvtss_f32(dir[1]));
        signs.z = copysign(1.0f, _mm_cvtss_f32(dir[2]));
    }
};
