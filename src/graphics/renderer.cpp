#include "renderer.h"

#include "graphics/lighting/sample.h"
#include "graphics/rays/frustum.h"
#include "graphics/primitives/basic/sphere.h"
#include "graphics/noise/gaussian.h"
#include "graphics/lighting/materials.h"
#include "graphics/noise/sampler.h"

#include "dev/gui.h"
#include "dev/debug.h"
#include "dev/profile.h"
#include <engine/physics/collision/sat.h>

void Renderer::init() {
    /* Try load the camera settings */
    FILE* f = fopen("camera.bin", "rb");
    if (f) {
        fread(&camera, 1, sizeof(Camera), f);
        fclose(f);
    }

#ifdef PROFILING
    const DWORD mask = 0b1 << (std::thread::hardware_concurrency() - 1);
    SetThreadAffinityMask(GetCurrentThread(), mask);

    profile_init(*this);
#endif

#ifdef DEV
    /* Assign the main camera */
    dev::main_camera = &camera;
#endif

#if DENOISE
    blur_in = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    blur_out = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
#endif

    /* Create the accumulator */
    accu = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (accu) memset(accu, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    albedo_buf = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (albedo_buf) memset(albedo_buf, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    prev_frame = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (prev_frame) memset(prev_frame, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));

    /* Load the blue noise sampler */
    bnoise = new BlueNoise();
}

/**
 * @brief Ray trace the scene.
 */
TraceResult Renderer::trace(Ray& ray, const u32 x, const u32 y, bool debug) const {
    /* Intersect the scene */
    const HitInfo hit = scene.intersect(ray);
    TraceResult result(hit);

    /* Handle special display modes, for debugging */
    if (dev::display_modes(hit, result)) return result;

    /* Return if we didn't hit a surface (result already has sky albedo) */
    if (hit.no_hit()) return result.no_reproject();

    /* Create a blue noise sampler */
    const NoiseSampler noise(bnoise, x, y, frame);

    /* Get the intersection point */
    const float3 hit_point = ray.intersection(hit);

    /* Evaluate material */
    MatEval eval;
    eval_material(eval, ray, hit, scene, noise);
    result.albedo = eval.albedo, result.irradiance = fmaxf(eval.irradiance, 0);

    return result;
}

TraceResult8x8 Renderer::trace(const RayPacket8x8& packet, const u32 x, const u32 y,
                               bool debug) const {
    TraceResult8x8 results;

    /* Intersect the scene */
    const PacketHit8x8 hits = scene.coherent_intersect(packet, debug);

    for (u32 r = 0; r < 8 * 8; r++) {
        const HitInfo& hit = hits.hits[r];
        TraceResult& result = results.results[r];
        result = TraceResult(hit);

        /* Handle special display modes, for debugging */
        if (dev::display_modes(hit, result)) continue;

        /* Return if we didn't hit a surface (result already has sky albedo) */
        if (hit.no_hit()) {
            result.no_reproject();
            continue;
        }

        /* Create a blue noise sampler */
        const NoiseSampler noise(bnoise, x + r % 8, y + r / 8, frame);

        /* Get the intersection point */
        const float3 hit_point = packet.rays[r].intersection(hit);

        /* Evaluate material */
        MatEval eval;
        eval_material(eval, packet.rays[r], hit, scene, noise);
        result.albedo = eval.albedo, result.irradiance = fmaxf(eval.irradiance, 0);
    }

    return results;
}

vector<float3> Renderer::path(const Ray& ray) const { 
    vector<float3> path;
    path.push_back(ray.origin);

    /* Intersect the scene */
    const HitInfo hit = scene.intersect(ray);

    /* Return if we didn't hit a surface (result already has sky albedo) */
    if (hit.no_hit()) {
        path.push_back(ray.origin + ray.dir * 1000.0f);
        return path;
    }

    path.push_back(ray.intersection(hit));

    Ray next_ray = ray;
    HitInfo next_hit = hit;
    for (u32 i = 0; i < 8; i++) {
        if (next_path_ray(next_ray, next_hit)) {
            /* Intersect the scene */
            next_hit = scene.intersect(next_ray);

            /* Return if we didn't hit a surface (result already has sky albedo) */
            if (next_hit.no_hit()) {
                path.push_back(next_ray.origin + next_ray.dir * 1000.0f);
                return path;
            }

            path.push_back(next_ray.intersection(next_hit));
        } else {
            return path;
        }
    }

    return path;
}

/**
 * @brief Called every frame.
 */
void Renderer::tick(Surface* screen, const f32 dt) {
    frame++;
    if (frame > 120) frame = 0;
    Timer t;

#if PACKET_TRACE
#ifndef PROFILING
#pragma omp parallel for schedule(dynamic)
#endif
    for (i32 y = 0; y < WIN_HEIGHT; y += 8) {
        for (i32 x = 0; x < WIN_WIDTH; x += 8) {
            const RayPacket8x8 packet = camera.get_packet8x8(x, y);
            const TraceResult8x8 results = trace(packet, x, y);

            for (i32 v = 0; v < 8; v++) {
                for (i32 u = 0; u < 8; u++) {
                    const i32 ix = x + u, iy = y + v;

                    const Ray& ray = packet.rays[v * 8 + u];
                    const TraceResult& r = results.results[v * 8 + u];
                    const float3 albedo = r.albedo;

                    /* Don't reproject if we hit the skybox */
                    if (r.depth >= BIG_F32 || not r.reproject) {
                        const float4 c = aces_approx(albedo);
                        screen->pixels[ix + iy * WIN_WIDTH] = RGBF32_to_RGB8(&c);
                        continue;
                    }

                    /* Accumulate and reproject */
                    const float3 ir = insert_accu(ix, iy, ray, r.irradiance, r.depth);
                    const float4 c = aces_approx(albedo * ir);
                    // albedo_buf[x + y * WIN_WIDTH] = albedo;
                    screen->pixels[ix + iy * WIN_WIDTH] = RGBF32_to_RGB8(&c);
                }
            }
        }
    }
#else
#ifndef PROFILING
#pragma omp parallel for schedule(dynamic)
#endif
    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            Ray ray = camera.get_primary_ray(x, y);

            /* Trace the scene */
            const TraceResult r = trace(ray, x, y);
            const float3 albedo = r.albedo;

            /* Don't reproject if we hit the skybox */
            if (r.depth >= BIG_F32 || not r.reproject) {
                const float4 c = aces_approx(albedo);
                screen->pixels[x + y * WIN_WIDTH] = RGBF32_to_RGB8(&c);
                continue;
            }

            /* Accumulate and reproject */
            const float3 ir = insert_accu(x, y, ray, r.irradiance, r.depth);
            const float4 c = aces_approx(albedo * ir);
            // albedo_buf[x + y * WIN_WIDTH] = albedo;
            screen->pixels[x + y * WIN_WIDTH] = RGBF32_to_RGB8(&c);
        }
    }
#endif

#if DENOISE
    /* Blur */
    memcpy(blur_in, accu, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    fast_gaussian_blur(blur_in, blur_out, WIN_WIDTH, WIN_HEIGHT, 2, 2.0f);

    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            u32 i = x + y * WIN_WIDTH;
            const float4 c = aces_approx(albedo_buf[x + y * WIN_WIDTH] * blur_out[i]);
            screen->pixels[i] = RGBF32_to_RGB8(&c);
        }
    }
#endif

    /* Swap the accumulator and the previous frame pointers */
    /* Faster than copying everything from the accu to the prev */
    float4* temp = accu;
    accu = prev_frame;
    prev_frame = temp;

#ifdef DEV
    dev::frame_time = t.elapsed();
#endif
}

void Renderer::gui(bool& running, const f32 dt) {
    // TODO: remove this
    // trace(dev::debug_packet, 0, 0, true);
    //db::draw_aabb(0, 1, 0xFFFF0000);
    // dev::debug_packet.setup_slice(0, 1, 32);
    // scene.cvv->intersect(dev::debug_packet, true);
    // dev::debug_packet.traverse(0, 1, 32);
}

void Renderer::shutdown() {
    /* Save the camera state */
    FILE* f = fopen("camera.bin", "wb");
    fwrite(&camera, 1, sizeof(Camera), f);
    fclose(f);

    delete bnoise;
}

/**
 * @brief Reproject onto the current frame and accumulate. (without tonemapping)
 * @return The color to display on screen for this pixel.
 */
inline float3 Renderer::insert_accu(const u32 x, const u32 y, const Ray& ray, const float3& c,
                                    const f32 d) const {
#ifdef DEV
    if (not dev::use_projection) return c;
#endif

    /* Reproject (c.w is the depth) */
    float3 acc_color = c;
    f32 confidence = 0.95f;

    const float2 prev_uv = camera.prev_pyramid.project(ray.origin + ray.dir * d);

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
        const f32 depth_diff = fabs(sample.w - (d + depth_delta));
        if (depth_diff < 0.1f) {
            confidence = max(confidence - depth_diff * 3.0f, 0.0f);
            acc_color = sample;
        }
    }

    /* Merge */
    const float3 color = (c * (1.0f - confidence)) + (acc_color * confidence);
    accu[x + y * WIN_WIDTH] = float4(color, d);
    return color;
}
