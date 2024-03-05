#pragma once

/* My own data types I like to use :) */
#include <stdint.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef float f32;
typedef double f64;

constexpr f32 BIG_F32 = 1e30f;
constexpr f32 EPSILON_F32 = std::numeric_limits<float>::epsilon();

/* SIMD data types */
#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>

typedef __m256 f256;
typedef __m256i i256;
typedef __m128 f128;
typedef __m128i i128;

const f128 BIG_F128 = _mm_set_ps1(BIG_F32);
