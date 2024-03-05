#pragma once

/* Basic SIMD friendly Ray structure. */
struct Ray {
    /* Ray origin */
    union {
        f128 o;
        struct {
            float3 origin;
            f32 _;
        };
    };
    /* Ray direction */
    union {
        f128 d;
        struct {
            float3 dir;
            f32 _;
        };
    };
    /* Ray direction reciprocal */
    union {
        f128 rd;
        struct {
            float3 r_dir;
            f32 _;
        };
    };
    /* Ray direction sign (-1 | 1) */
    union {
        f128 sd;
        struct {
            float3 sign_dir;
            f32 _;
        };
    };
    f32 t = BIG_F32;
    u32 steps = 0;

    Ray() = default;
    inline Ray(const float3& origin, const float3& dir)
        : origin(origin), dir(dir), r_dir(1.0f / dir), sign_dir(sign_of_dir(dir)) {}

    inline f32 intersect_aabb(const f128& min, const f128& max) const {
        /* Idea to use fmsub to save 1 instruction came from
         * <http://www.joshbarczak.com/blog/?p=787> */
        const f128 ord = _mm_mul_ps(o, rd);
        const f128 t1 = _mm_fmsub_ps(min, rd, ord);
        const f128 t2 = _mm_fmsub_ps(max, rd, ord);

        /* Find the near and far intersection point */
        const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);

        /* Set the 4th element to 0 */
        const f128 tmin4 = (f128&)_mm_slli_si128((i128&)vmin4, 4);

        /* Get the horizontal minimum and maximum "t" */
        const f32 tmax = _mm_hmin3_ps(vmax4), tmin = _mm_hmax_ps(tmin4);

        f32 r = BIG_F32;
        if (tmax >= tmin) {
            r = tmin;
        }
        return r;
    }

   private:
    __forceinline i32 getsign(const f32 f) const { return (i32)(((u32&)f) >> 31) * 2 - 1; }

    /* Get the signs of a vector (-1 | 1) */
    __forceinline float3 sign_of_dir(float3 d) const {
        {
            float3 signs = float3();
#if 1
            signs.x = -getsign(dir.x);
            signs.y = -getsign(dir.y);
            signs.z = -getsign(dir.z);
#else
            signs.x = copysign(1.0f, dir.x);
            signs.y = copysign(1.0f, dir.y);
            signs.z = copysign(1.0f, dir.z);
#endif
            return signs;
        }
    }
};
