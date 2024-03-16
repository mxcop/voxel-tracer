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

constexpr u32 BVH_SHAPES = 9;

struct TraceResult {
    float4 albedo = 0, irradiance = 0;
};

class Renderer : public TheApp {
   public:
    /** @brief Initialize the application */
    void init();
    /** @brief Trace the scene */
    TraceResult trace(Ray& ray, HitInfo& hit, const u32 x, const u32 y) const;
    /** @brief Called as often as possible */
    void tick(f32 dt);
    /** @brief Called after "tick" for ImGUI rendering */
    void gui(f32 dt);
    /** @brief Called before the application closes */
    void shutdown();

    /* User input */
    void MouseUp(int button) {}
    void MouseDown(int button);
    void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
    void MouseWheel(float y) {}
    void KeyUp(int key) {}
    void KeyDown(int key) {}

    /**
     * @brief Reproject onto the current frame and accumulate. (with tonemapping)
     * @return The color to display on screen for this pixel.
     */
    inline float4 insert_accu(const u32 x, const u32 y, const Ray& ray, const float4& c) const;

    /**
     * @brief Reproject onto the current frame and accumulate. (without tonemapping)
     * @return The color to display on screen for this pixel.
     */
    inline float4 insert_accu_raw(const u32 x, const u32 y, const Ray& ray, const float4& c) const;

    int2 mousePos;
    u32 frame = 0u;
    float3 sun_dir = {-0.619501f, 0.465931f, -0.631765f};

    /* Game camera */
    Camera camera;
    f32 depth_delta = 0;

    /* Lights */
    vector<LightSource> lights;
    vector<SphereLight> area_lights;

    /* Scene */
#if USE_BVH
    Traceable* shapes[BVH_SHAPES];
    Bvh* bvh = nullptr;
    Sphere* test_vv = nullptr;
    Sphere* test_plane_vv = nullptr;
    OVoxelVolume* arm_vv[4] = {};
#else
    VoxelVolume* volume = nullptr;
    // BrickVolume* volume = nullptr;
#endif

    /* Textures */
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

    /* Physics testing */
    PhyWorld world;
    PhyObject* test_obj = nullptr;
    PhyObject* test_plane_obj = nullptr;
    RobotArm arm;
};
