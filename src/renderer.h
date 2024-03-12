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

class Renderer : public TheApp {
   public:
    void init();
    u32 trace(Ray& ray, const u32 x, const u32 y) const;
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
    inline float4 insert_accu(const u32 x, const u32 y, const float4& c) const {
        return aces_approx(insert_accu_raw(x, y, c));
    }

    /* Insert a new color into the accumulator without tonemapping (returns the color to display) */
    inline float4 insert_accu_raw(const u32 x, const u32 y, const float4& c) const {
        const float4 new_color = c;
        const float4 acc_color = accu[x + y * WIN_WIDTH];
        if (accu_reset) {
            accu[x + y * WIN_WIDTH] = new_color;
            return new_color;
        }
        const float4 color = new_color * (0.1f) + acc_color * (0.9f);
        accu[x + y * WIN_WIDTH] = color;
        return color;
    }

    /* Reset the accumulator */
    inline void reset_accu() {
        accu_len = 1u, memset(accu, 0, (size_t)WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
        frame = 0;
        accu_reset = true;
    };

    int2 mousePos;
    f32 frame_time = 1.0f;
    u32 frame = 0u;
    float3 sun_dir = {-0.619501f, 0.465931f, -0.631765f};

    Camera camera;
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

    /* Accumulator */
    float4* accu = nullptr;
    mutable u32 accu_len = 1u;
    mutable bool accu_reset = true;

    bool fast_mode = true;

    /* Physics testing */
    PhyWorld world;
    PhyObject* test_obj = nullptr;
    PhyObject* test_plane_obj = nullptr;
};
