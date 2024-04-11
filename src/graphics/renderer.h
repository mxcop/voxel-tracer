#pragma once

#include "graphics/camera.h"
#include "graphics/noise/blue.h"
#include "graphics/skydome.h"
#include "graphics/lighting/sphere-light.h"
#include "graphics/bvh.h"
#include "graphics/primitives/vv.h"
#include "graphics/tonemap.h"
#include "graphics/scene.h"
#include "graphics/lighting/materials.h"

#include "engine/physics/world.h"

/* Quite expensive! */
#define DENOISE 0

struct TraceResult {
    float4 albedo = 0, irradiance = 0;
    f32 depth = BIG_F32;
    bool reproject = true;

    TraceResult() = default;
    explicit TraceResult(const HitInfo& hit)
        : albedo(hit.albedo), irradiance(0), depth(hit.depth), reproject(true) {}

    /** @brief Don't reproject this sample (returns itself) */
    TraceResult& no_reproject() {
        reproject = false;
        return *this;
    };
};

struct TraceResult8x8 {
    TraceResult results[8 * 8];
};

class Renderer {
   public:
    /** @brief Initialize the application */
    void init();
    /** @brief Trace the scene */
    TraceResult trace(Ray& ray, const u32 x, const u32 y) const;
    TraceResult8x8 trace(const RayPacket8x8& packet, const u32 x, const u32 y, bool debug = false) const;
    vector<float3> path(const Ray& ray) const;
    /** @brief Called as often as possible */
    void tick(Surface* screen, const f32 dt);
    /** @brief Called after "tick" for ImGUI rendering */
    void gui(bool& running, const f32 dt);
    /** @brief Called before the application closes */
    void shutdown();

    /**
     * @brief Reproject onto the current frame and accumulate. (without tonemapping)
     * @return The color to display on screen for this pixel.
     */
    inline float3 insert_accu(const u32 x, const u32 y, const Ray& ray, const float3& c,
                                  const f32 d) const;

    u32 frame = 0u;

    /* Game camera */
    Camera camera;
    f32 depth_delta = 0;

    /* Scene */
    Scene scene;

    /* Blue noise sampler */
    const BlueNoise* bnoise = nullptr;

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
    //PhyWorld world;
    //PhyObject* test_obj = nullptr;
    //PhyObject* test_plane_obj = nullptr;
};
