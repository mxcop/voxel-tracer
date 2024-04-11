// Template, IGAD version 3
// IGAD/NHTV/UU - Jacco Bikker - 2006-2022
#pragma once

/* Settings */
#if 1
constexpr u32 WIN_WIDTH = 1280;
constexpr u32 WIN_HEIGHT = 720;
#elif 0
constexpr u32 WIN_WIDTH = 1920;
constexpr u32 WIN_HEIGHT = 1080;
#else
constexpr u32 WIN_WIDTH = 256;
constexpr u32 WIN_HEIGHT = 212;
#endif
constexpr float2 WIN_SIZE = float2(WIN_WIDTH, WIN_HEIGHT);

constexpr f32 VOXEL = 1.0f / 20.0f;

/* Switch between a BVH world and a voxel volume */
#define USE_BVH 1

/* For profiling purposes */
//#define PROFILING

/* Use packet tracing for primary rays */
#define PACKET_TRACE 0

/* Use more accurate pyramid tracing */
#define ACCURATE_PYRAMID_TRACING 1

/* Inline functions */
static inline f32 _min(const f32 a, const f32 b) { return a < b ? a : b; };
static inline f32 _max(const f32 a, const f32 b) { return a > b ? a : b; };

/* Some constants */
#define PI 3.14159265358979323846264f
#define INVPI 0.31830988618379067153777f
#define INV2PI 0.15915494309189533576888f
#define TWOPI 6.28318530717958647692528f
#define FOURPI 12.5663706144f
#define SQRT_PI_INV 0.56418958355f
#define LARGE_FLOAT 1e34f

// IMPORTANT NOTE ON OPENCL COMPATIBILITY ON OLDER LAPTOPS:
// Without a GPU, a laptop needs at least a 'Broadwell' Intel CPU (5th gen, 2015):
// Intel's OpenCL implementation 'NEO' is not available on older devices.
// Same is true for Vulkan, OpenGL 4.0 and beyond, as well as DX11 and DX12.