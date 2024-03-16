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
#include "graphics/primitives/basic/sphere.h"

#include "engine/physics/world.h"
#include "game/robot-arm.h"

/* Quite expensive! */
#define DENOISE 0

/* What happened?? */
#define ANGRY_MODE 0

constexpr u32 BVH_SHAPES = 9;

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

            constexpr f32 MAX_U = 1.0f - (1.0f / WIN_WIDTH) * 2;
            constexpr f32 MAX_V = 1.0f - (1.0f / WIN_HEIGHT) * 2;

            if (prev_uv.x > 0 && prev_uv.x < MAX_U && prev_uv.y > 0 && prev_uv.y < MAX_V) {
                /* Bilinear Sampling */
                const float2 center = prev_uv * WIN_SIZE + 0.5f;
                
                /* Sample points */
                const float2 tl_p = prev_uv * WIN_SIZE;
                const float2 tr_p = prev_uv * WIN_SIZE + float2(1, 0);
                const float2 bl_p = prev_uv * WIN_SIZE + float2(0, 1);
                const float2 br_p = prev_uv * WIN_SIZE + float2(1, 1);

                /* Center pixel */
                const float2 center_p = floorf(center + 0.5f);

                /* Sample weights */
                const f32 tl_w = fabs((tl_p.x - center_p.x) * (tl_p.y - center_p.y));
                const f32 tr_w = fabs((tr_p.x - center_p.x) * (tr_p.y - center_p.y));
                const f32 bl_w = fabs((bl_p.x - center_p.x) * (bl_p.y - center_p.y));
                const f32 br_w = 1.0f - (tl_w + tr_w + bl_w);

                /* Fetch the samples */
                const float3 tl_s = prev_frame[(i32)tl_p.x + (i32)tl_p.y * WIN_WIDTH] * tl_w;
                const float3 tr_s = prev_frame[(i32)tr_p.x + (i32)tr_p.y * WIN_WIDTH] * tr_w;
                const float3 bl_s = prev_frame[(i32)bl_p.x + (i32)bl_p.y * WIN_WIDTH] * bl_w;
                const float3 br_s = prev_frame[(i32)br_p.x + (i32)br_p.y * WIN_WIDTH] * br_w;

                /* Merge the samples, and use the center depth value */
                const f32 depth = prev_frame[(i32)center.x + (i32)center.y * WIN_WIDTH].w;
                const float4 sample = float4(tl_s + tr_s + bl_s + br_s, depth);

                /* Depth rejection (take into account camera movement "depth_delta") */
                const f32 depth_diff = fabs(sample.w - (c.w + depth_delta));
                if (depth_diff < 0.2f) {
                    confidence = max(confidence - depth_diff * 3.0f, 0.0f);
                    acc_color = sample;
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
    Traceable* shapes[BVH_SHAPES];
    Bvh* bvh = nullptr;
    Sphere* test_vv = nullptr;
    Sphere* test_plane_vv = nullptr;
    // OVoxelVolume* arm_vv = nullptr;
    OVoxelVolume* arm_vv[4] = {};
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
    RobotArm arm;
};
