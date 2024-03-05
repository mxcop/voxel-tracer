#pragma once

/* Holds information about a ray hit. */
struct HitInfo {
    f32 depth = BIG_F32;
    float4 albedo = {};
    float3 normal = {};
    u32 steps = 0; /* For debugging! */
};

struct PacketHitInfo {
    f128 depth;
    f128 exit_t;
    u32 steps = 0; /* For debugging! */
};
