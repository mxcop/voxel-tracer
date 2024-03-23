#pragma once

//struct ext {
//    f32 min, max;
//};

/* Coherent ray packet 4x4, cache line aligned. */
struct alignas(64) CoherentPacket8x8 {
    /* [du_min, du_max, dv_min, dv_max] */
    /* Factor by which the slice grows each step. */
    union {
        f128 delta_slice;
        struct {
            f32 du_min, du_max, dv_min, dv_max;
        };
    };
    /* [u_min, u_max, v_min, v_max] */
    /* Extend of the current slice. */
    union {
        f128 slice;
        struct {
            f32 u_min, u_max, v_min, v_max;
        };
    };
    /* Traversal axis. */
    u32 k, u, v;
    /* Distance along major axis. */
    f32 k_t = 0;

    /* Packet origin point. */
    float3 origin;
    /* Packet ray directions. */
    float3 rays[4 * 4];
    /* Packet direction signs. */
    // float3 signs;

    void setup_slice(const float3& min, const float3& max, const f32 vpu);

    void traverse(const float3& min, const float3& max, const f32 vpu);

   private:
    // FOR TESTING ONLY
    void draw_slice(const f32 vpu) const;

    __forceinline i32 getsign(const f32 f) const { return (i32)(((u32&)f) >> 31) * 2 - 1; }
    f32 entry(const f32 ro, const f32 rd, const f32 min, const f32 max) const;
};

/* SIMD (SSE) Ray packet structure. */
struct alignas(64) RayPacket128 {
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

    RayPacket128() = default;
    inline RayPacket128(const f128 shared_origin[3], const f128 dir[3]) {
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

    /* Fast ray packet vs AABB intersection test, using SSE, FMA */
    inline f128 intersect_aabb(const float3& min, const float3& max) const {
        f128 tmin = _mm_setzero_ps(), tmax = BIG_F128;

        for (u32 a = 0; a < 3; ++a) {
            const f128 bmin = _mm_set_ps1(min[a]);
            const f128 bmax = _mm_set_ps1(max[a]);

            const f128 ord = _mm_mul_ps(ro[a], ird[a]);
            const f128 t1 = _mm_fmsub_ps(bmin, ird[a], ord);
            const f128 t2 = _mm_fmsub_ps(bmax, ird[a], ord);

            tmin = _mm_min_ps(_mm_max_ps(t1, tmin), _mm_max_ps(t2, tmin));
            tmax = _mm_max_ps(_mm_min_ps(t1, tmax), _mm_min_ps(t2, tmax));
        }

        /* Use a mask to remove non-intersections (tmin > tmax) */
        const f128 mask = _mm_cmple_ps(tmin, tmax);
        return _mm_blendv_ps(BIG_F128, tmin, mask);
    }
};

/* SIMD (AVX) Ray packet structure. */
struct alignas(64) RayPacket256 {
    union {
        f256 ro[3]; /* Ray origin (XYZ) */
        struct {
            f256 ro_x;
            f256 ro_y;
            f256 ro_z;
        };
    };
    union {
        f256 rd[3]; /* Ray direction (XYZ) */
        struct {
            f256 rd_x;
            f256 rd_y;
            f256 rd_z;
        };
    };
    union {
        f256 ird[3]; /* Ray direction reciprocal (XYZ) */
        struct {
            f256 ird_x;
            f256 ird_y;
            f256 ird_z;
        };
    };
    float3 signs; /* Ray signs are shared by all rays in packet! */

    RayPacket256() = default;
    inline RayPacket256(const f256 shared_origin[3], const f256 dir[3]) {
        ro[0] = shared_origin[0];
        ro[1] = shared_origin[1];
        ro[2] = shared_origin[2];
        rd[0] = dir[0];
        rd[1] = dir[1];
        rd[2] = dir[2];
        ird[0] = _mm256_rcp_ps(dir[0]);
        ird[1] = _mm256_rcp_ps(dir[1]);
        ird[2] = _mm256_rcp_ps(dir[2]);
        signs.x = copysign(1.0f, _mm256_cvtss_f32(dir[0]));
        signs.y = copysign(1.0f, _mm256_cvtss_f32(dir[1]));
        signs.z = copysign(1.0f, _mm256_cvtss_f32(dir[2]));
    }
};
