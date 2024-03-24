#pragma once

/* Holds information about a ray hit. */
struct HitInfo {
    f32 depth = BIG_F32;
    u8 material = 0;
    float4 albedo = {};
    float3 normal = {};
    u32 steps = 0; /* For debugging! */

    /** @return True if no hit occured. */
    bool no_hit() const { return depth == BIG_F32; };
};

struct PacketHitInfo {
    f128 depth;
    f128 exit_t;
    u32 steps = 0; /* For debugging! */
};

struct alignas(64) CoherentHit4x4 {
    f32 depth[4 * 4] = {};

    CoherentHit4x4() {
        for (u32 r = 0; r < 4*4; r++) {
            depth[r] = BIG_F32;
        }
    }
};
