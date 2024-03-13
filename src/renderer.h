#pragma once

#include "graphics/camera.h"
#include "graphics/primitives/voxel-volume.h"
#include "graphics/primitives/brick-volume.h"
#include "graphics/light.h"
#include "graphics/noise/blue.h"
#include "graphics/skydome.h"
#include "graphics/lighting/sphere-light.h"
#include "graphics/bvh.h"
#include "graphics/primitives/vv.h"
#include "graphics/tonemap.h"

#include "engine/physics/world.h"
#include <graphics/primitives/basic/sphere.h>

/* Quite expensive! */
#define DENOISE 0

struct TraceResult {
    float4 albedo = 0, irradiance = 0;
};

class Renderer : public TheApp {
   public:
    void init();
    TraceResult trace(Ray& ray, HitInfo& hit, const u32 x, const u32 y) const;
    void tick(f32 dt);
    void gui(f32 dt);
    void shutdown();
    /* User input */
    void MouseUp(int button) {}
    void MouseDown(int button);
    void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
    void MouseWheel(float y) {}
    void KeyUp(int key) {}
    void KeyDown(int key) {}

    /* Insert a new color into the accumulator (returns the color to display) */
    inline float4 insert_accu(const u32 x, const u32 y, const Ray& ray, const float4& c) const {
        return aces_approx(insert_accu_raw(x, y, ray, c));
    }

    /* Insert a new color into the accumulator without tonemapping (returns the color to display) */
    inline float4 insert_accu_raw(const u32 x, const u32 y, const Ray& ray, const float4& c) const {
        /* Reproject (color.w is the depth) */
        float4 acc_color = c;
        f32 confidence = 0.9f;
        if (c.w < BIG_F32) {
            const float2 prev_uv = camera.prev_pyramid.project(ray.origin + ray.dir * c.w);

            constexpr f32 MAX_U = 1.0f - (1.0f / WIN_WIDTH) * 0.5f;
            constexpr f32 MAX_V = 1.0f - (1.0f / WIN_HEIGHT) * 0.5f;

            if (prev_uv.x >= 0 && prev_uv.x < MAX_U && prev_uv.y >= 0 && prev_uv.y < MAX_V) {
                const i32 px = static_cast<i32>((prev_uv.x * WIN_WIDTH) + 0.5f);
                const i32 py = static_cast<i32>((prev_uv.y * WIN_HEIGHT) + 0.5f);

                const float4 prev_c = prev_frame[px + py * WIN_WIDTH];

                /* Depth rejection (take into account camera movement "depth_delta") */
                const f32 depth_diff = fabs(prev_c.w - (c.w + depth_delta));
                if (depth_diff < 0.2f) {
                    confidence = max(confidence - depth_diff * 3.0f, 0.0f);
                    acc_color = prev_c;
                }
            }
        }

        /* Merge */
        const float4 color = (c * (1.0f - confidence)) + (acc_color * confidence);
        accu[x + y * WIN_WIDTH] = float4(color.x, color.y, color.z, c.w);
        return color;
    }

    int2 mousePos;
    f32 frame_time = 1.0f;
    u32 frame = 0u;
    float3 sun_dir = {-0.619501f, 0.465931f, -0.631765f};

    Camera camera;
    f32 depth_delta = 0;
    vector<LightSource> lights;
    vector<SphereLight> area_lights;
    Surface* texture = nullptr;

#if USE_BVH
    Traceable* shapes[7];
    Bvh* bvh = nullptr;
    Sphere* test_vv = nullptr;
    Sphere* test_plane_vv = nullptr;
    OVoxelVolume* arm_vv = nullptr;
#else
    VoxelVolume* volume = nullptr;
    // BrickVolume* volume = nullptr;
#endif

    BlueNoise* bnoise = nullptr;
    SkyDome skydome;

    /* === BUFFERS === */
    /* Accumulator (reprojection) */
    float4* accu = nullptr;
    float4* albedo_buf = nullptr;
    /* Stores XYZ: color, W: depth */
    float4* prev_frame = nullptr;

#if DENOISE
    float4 *blur_in = nullptr, *blur_out = nullptr;
#endif

    bool fast_mode = true;

    /* Physics testing */
    PhyWorld world;
    PhyObject* test_obj = nullptr;
    PhyObject* test_plane_obj = nullptr;
};
