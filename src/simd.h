#pragma once

/* Contains 4, 3 component vectors, stored in vertical fashion. */
union f128q {
    f128 q[3] = {};
    struct {
        f128 x, y, z;
    };
};

/* Contains 4, 3 component vectors, stored in vertical fashion. */
union i128q {
    i128 q[3] = {};
    struct {
        i128 x, y, z;
    };
};

/* SSE, Clamp and then floor. */
inline static i128 _mm_floorclamp_ps(const f128 m, const i128 min, const i128 max) {
    return _mm_min_epi32(_mm_max_epi32(_mm_cvtps_epi32(_mm_floor_ps(m)), min), max);
}

__forceinline static i128 _mm_clamp_epi32(const i128 v, const u32 min, const u32 max) {
    return _mm_min_epi32(_mm_max_epi32(v, _mm_set1_epi32(min)), _mm_set1_epi32(max));
}

/* SSE, Absolute float vector. */
inline static f128 _mm_abs_ps(const f128 m) { return _mm_andnot_ps(_mm_set_ps1(-0.0f), m); }

/* SSE2 128 bit horizontal maximum */
__forceinline f32 _mm_hmax_ps(const f128 v) {
    const i128 x = _mm_shuffle_epi32((i128&)v, 0b01'00'11'10);
    const f128 y = _mm_max_ps(v, (f128&)x); /* 4 -> 2 */
    const i128 z = _mm_shuffle_epi32((i128&)y, 0b00'01'00'01);
    const f128 r = _mm_max_ps(y, (f128&)z); /* 2 -> 1 */
    return _mm_cvtss_f32(r);
}

/* SSE2 128 bit horizontal minimum (ignores 4th element) */
__forceinline f32 _mm_hmin3_ps(const f128 v) {
    i128 x = _mm_shuffle_epi32((i128&)v, 0b01'00'00'10);
    f128 y = _mm_min_ps(v, (f128&)x); /* 4 -> 2 */
    i128 z = _mm_shuffle_epi32((i128&)y, 0b00'01'00'01);
    f128 r = _mm_min_ps(y, (f128&)z); /* 2 -> 1 */
    return _mm_cvtss_f32(r);
}

/* SSE2 128 bit horizontal minimum */
__forceinline f32 _mm_hmin_ps(const f128 v) {
    const i128 x = _mm_shuffle_epi32((i128&)v, 0b01'00'11'10);
    const f128 y = _mm_min_ps(v, (f128&)x); /* 4 -> 2 */
    const i128 z = _mm_shuffle_epi32((i128&)y, 0b00'01'00'01);
    const f128 r = _mm_min_ps(y, (f128&)z); /* 2 -> 1 */
    return _mm_cvtss_f32(r);
}
