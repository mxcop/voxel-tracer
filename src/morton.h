#pragma once

/* If on an x64 processor, only works with BMI2 */
#if defined(_M_X64) && (defined(__BMI2__) || defined(__AVX2__))

/* 3D Morton BMI masks */
constexpr u64 BMI_3D_X_MASK = 0x9249249249249249;
constexpr u64 BMI_3D_Y_MASK = 0x2492492492492492;
constexpr u64 BMI_3D_Z_MASK = 0x4924924924924924;
constexpr u64 BMI_3D_MASKS[3] = {BMI_3D_X_MASK, BMI_3D_Y_MASK, BMI_3D_Z_MASK};

/* 3D coordinate to morton code. */
inline u64 morton_encode(const u32 x, const u32 y, const u32 z) {
    return _pdep_u64((u64)x, BMI_3D_X_MASK) | _pdep_u64((u64)y, BMI_3D_Y_MASK) |
           _pdep_u64((u64)z, BMI_3D_Z_MASK);
}

/* Fast version of morton encode for SSE using BMI2 (limited to 11/10 bits per axis) */
inline i128 morton_encode(const i128 x, const i128 y, const i128 z) {
    /* Bit deposit masks */
    constexpr u64 X_MASK = 0x4924924949249249; /* 11 bits */
    constexpr u64 Y_MASK = 0x9249249292492492; /* 11 bits */
    constexpr u64 Z_MASK = 0x2492492424924924; /* 10 bits */

    /* Adjust bit format for conversion */
    const i128 sx = _mm_srli_epi64(x, (32 - 11 /* bits */));
    const i128 sy = _mm_srli_epi64(y, (32 - 11 /* bits */));
    const i128 sz = _mm_srli_epi64(z, (32 - 10 /* bits */));
    const i128 b11m = _mm_set1_epi64x(0b11111111111);
    const i128 b10m = _mm_set1_epi64x(0b01111111111);
    /* TIP: technically I could 'not' do the 'andnot' it will only fail with very large numbers */
    const i128 cx = _mm_or_si128(_mm_andnot_si128(b11m, sx), _mm_and_si128(b11m, x));
    const i128 cy = _mm_or_si128(_mm_andnot_si128(b11m, sy), _mm_and_si128(b11m, y));
    const i128 cz = _mm_or_si128(_mm_andnot_si128(b10m, sz), _mm_and_si128(b10m, z));

    /* Gather lower and upper 64 bit integers */
    const u64 x1 = _mm_cvtsi128_si64(cx);
    const u64 y1 = _mm_cvtsi128_si64(cy);
    const u64 z1 = _mm_cvtsi128_si64(cz);
    const u64 x2 = _mm_extract_epi64(cx, 1);
    const u64 y2 = _mm_extract_epi64(cy, 1);
    const u64 z2 = _mm_extract_epi64(cz, 1);

    /* Perform parallel bit deposits */
    const u64 p1 = _pdep_u64(x1, X_MASK);
    const u64 p2 = _pdep_u64(x2, X_MASK);
    const u64 p3 = _pdep_u64(y1, Y_MASK);
    const u64 p4 = _pdep_u64(y2, Y_MASK);
    const u64 p5 = _pdep_u64(z1, Z_MASK);
    const u64 p6 = _pdep_u64(z2, Z_MASK);

    // x | y | z
    /* Combine the results */
    const u64 m135 = p1 | p3 | p5;
    const u64 m246 = p2 | p4 | p6;

    /* Store and return the resulting vector */
    return _mm_set_epi64x(m246, m135);
}

inline i128 morton_encode_ex(const i128 x, const i128 y, const i128 z) {
    return _mm_set_epi32(morton_encode(x.m128i_u32[3], y.m128i_u32[3], z.m128i_u32[3]),
                         morton_encode(x.m128i_u32[2], y.m128i_u32[2], z.m128i_u32[2]),
                         morton_encode(x.m128i_u32[1], y.m128i_u32[1], z.m128i_u32[1]),
                         morton_encode(x.m128i_u32[0], y.m128i_u32[0], z.m128i_u32[0]));
}

/* Morton code to 3D coordinate. */
inline void morton_decode(const u64 m, u32& x, u32& y, u32& z) {
    x = (u32)(_pext_u64(m, BMI_3D_X_MASK));
    y = (u32)(_pext_u64(m, BMI_3D_Y_MASK));
    z = (u32)(_pext_u64(m, BMI_3D_Z_MASK));
}

/* 2D Morton BMI masks */
constexpr u64 BMI_2D_X_MASK = 0x5555555555555555;
constexpr u64 BMI_2D_Y_MASK = 0xAAAAAAAAAAAAAAAA;

/* 2D coordinate to morton code. */
inline u64 morton_encode(const u32 x, const u32 y) {
    return _pdep_u64((u64)x, BMI_2D_X_MASK) | _pdep_u64((u64)y, BMI_2D_Y_MASK);
}

/* Morton code to 2D coordinate. */
inline void morton_decode(const u64 m, u32& x, u32& y) {
    x = (u32)(_pext_u64(m, BMI_2D_X_MASK));
    y = (u32)(_pext_u64(m, BMI_2D_Y_MASK));
}

#else

/* Seperate bits from a given integer 3 positions apart. */
inline u64 split3(u32 a) {
    /* We only care about the first 21 bits */
    u64 x = a & 0x1fffff;
    /* 0000000000011111000000000000000000000000000000001111111111111111 */
    x = (x | x << 32) & 0x1f00000000ffff;
    /* 0000000000011111000000000000000011111111000000000000000011111111 */
    x = (x | x << 16) & 0x1f0000ff0000ff;
    /* 0001000000001111000000001111000000001111000000001111000000000000 */
    x = (x | x << 8) & 0x100f00f00f00f00f;
    /* 0001000011000011000011000011000011000011000011000011000100000000 */
    x = (x | x << 4) & 0x10c30c30c30c30c3;
    /* 0001001001001001001001001001001001001001001001001001001001001001 */
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

/* 3D coordinate to morton code. */
inline u64 morton_encode(const u32 x, const u32 y, const u32 z) {
    return split3(x) << 0 | split3(y) << 1 | split3(z) << 2;
}

/* Basically does the inverse of what "split3" does. */
inline u32 third_bits(const u64 m) {
    /* 0001001001001001001001001001001001001001001001001001001001001001 */
    u64 x = m & 0x1249249249249249;
    /* 0001000011000011000011000011000011000011000011000011000100000000 */
    x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3;
    /* 0001000000001111000000001111000000001111000000001111000000000000 */
    x = (x ^ (x >> 4)) & 0x100f00f00f00f00f;
    /* 0000000000011111000000000000000011111111000000000000000011111111 */
    x = (x ^ (x >> 8)) & 0x1f0000ff0000ff;
    /* 0000000000011111000000000000000000000000000000001111111111111111 */
    x = (x ^ (x >> 16)) & 0x1f00000000ffff;
    /* We only care about the first 21 bits */
    x = (x ^ (x >> 32)) & 0x1fffff;
    return (u32)x;
}

/* Morton code to 3D coordinate. */
inline void morton_decode(const u64 m, u32& x, u32& y, u32& z) {
    x = third_bits(m), y = third_bits(m >> 1), z = third_bits(m >> 2);
}

#endif
