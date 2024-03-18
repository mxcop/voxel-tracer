#include "graphics/lighting/sample.h"
#include "graphics/rays/frustum.h"
#include "graphics/primitives/basic/sphere.h"
#include <graphics/noise/gaussian.h>

#include "dev/gui.h"
#include "dev/debug.h"
#include "dev/profile.h"

u32 HSBtoRGB(f32 h, f32 s, f32 v) {
    assert(-360 <= h && h <= 360 && "h must be within [-360; 360]");
    assert(0 <= s && s <= 1 && "s must be within [0; 100]");
    assert(0 <= v && v <= 1 && "v must be within [0; 100]");

    // Keep h within [0; 359], s and v within [0; 1]
    if (h >= 360) h -= 360;
    if (h < 0) h += 360;

    // Convert hsv to rgb. This algorithm is described at
    // https://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV
    const f32 C = v * s;
    const f32 hd = h / 60;
    const int hi = int(hd);
    const f32 hd_mod2 = hd - hi + hi % 2;
    const f32 X = C * (1 - fabs(hd_mod2 - 1));
    f32 r, g, b;

    switch (hi) {
        case 0:
            r = C;
            g = X;
            b = 0;
            break;
        case 1:
            r = X;
            g = C;
            b = 0;
            break;
        case 2:
            r = 0;
            g = C;
            b = X;
            break;
        case 3:
            r = 0;
            g = X;
            b = C;
            break;
        case 4:
            r = X;
            g = 0;
            b = C;
            break;
        case 5:
            r = C;
            g = 0;
            b = X;
            break;
        // This should never happen
        default:
            return 0;
    }

    const f32 m = v - C;
    r += m;
    g += m;
    b += m;

    // Scale r, g, b to [0; 255] and pack them into a number 0xRRGGBB
    return (int(r * 255) << 16) + (int(g * 255) << 8) + int(b * 255);
}

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
#endif

#ifdef DEV
    /* Assign the main camera */
    dev::main_camera = &camera;
#endif

    /* Create the accumulator */
    accu = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (accu) memset(accu, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    albedo_buf = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (albedo_buf) memset(albedo_buf, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    prev_frame = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (prev_frame) memset(prev_frame, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));

#if DENOISE
    blur_in = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    blur_out = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
#endif

    bnoise = new BlueNoise();
    skydome = SkyDome("assets/kiara_1_dawn_4k.hdr");

#if USE_BVH
    const f32 VOXEL = 1.0f / 20.0f;

    /* Crates */
    shapes[0] = new OVoxelVolume({0, 0, 0}, "assets/vox/material-test.vox");
    //shapes[1] =
    //    new OVoxelVolume(float3(-3.0f, 2.5f - VOXEL * 3, -0.5f), "assets/vox/crate-16h.vox");
    //shapes[2] = new OVoxelVolume(float3(1.0f + VOXEL * 17, 2.5f, -0.5f), "assets/vox/crate-16.vox");

    bvh = new Bvh(BVH_SHAPES, shapes);

    /* Physics testing */
    // test_obj = world.add_object(PhyObject(new SphereCollider(0, 1.5f), float3(0, 10, 1.0f),
    // 0.01f)); test_plane_obj =
    //     world.add_object(PhyObject(new BoxCollider(float3(-4, -1, -4), float3(4, 0, 4)), 0, 0));
    
#else
    volume = new VoxelVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
#endif

#ifdef PROFILING
    profile_init(*this);
#endif
}

TraceResult Renderer::trace(Ray& ray, HitInfo& hit, const u32 x, const u32 y) const {
    TraceResult result;
#if USE_BVH
    const Bvh* volume = this->bvh;
#endif
    hit = volume->intersect(ray);

    float3 hit_pos = ray.origin + ray.dir * hit.depth + hit.normal * 0.00001f;

    /* Reflections */
#if 1
    constexpr u32 MAX_BOUNCES = 8;
    /* TODO: also keep looping if the material hit is transmissive */
    for (u32 i = 0; i < MAX_BOUNCES && hit.albedo.w > 0.0f; ++i) {
        /* TODO: distingus between glass and mirror! */
        /* Blue noise + R2 (cosine weighted distribution) */
        float3 raw_noise = bnoise->sample_3d(x, y);
        float3 quasi_noise;
        quasi_noise.x = fmod(raw_noise.x + R2X * (f32)frame, 1.0f);
        quasi_noise.y = fmod(raw_noise.y + R2Y * (f32)frame, 1.0f);
        quasi_noise.z = fmod(raw_noise.z + R2Z * (f32)frame, 1.0f);
        const float3 jitter = (quasi_noise * hit.albedo.w - (hit.albedo.w * 0.5f));
        const float3 reflect_dir = normalize(reflect(ray.dir, hit.normal) + jitter);

        ray = Ray(hit_pos, reflect_dir);
        const u32 steps = hit.steps;
        hit = volume->intersect(ray);
        hit.steps += steps; /* <- carry step count */

        /* Break if the ray missed */
        if (hit.depth == BIG_F32) {
            break;
        }

        /* Update the intersection point */
        hit_pos = ray.origin + ray.dir * hit.depth /* + hit.normal * 0.000001f*/;

        if (hit.material > 0 && hit.material <= 8) {
            ray = Ray(hit_pos, ray.dir);
            ray.medium_id = hit.material;
            hit = volume->intersect(ray);
        }
    }
#endif
    result.albedo = float4(0);

#if 0
    /* Textures */
    f32 u, v;

    const float3 p = hit_pos / 8.0f;
    if (hit.normal.x) {
        u = p.y - floorf(p.y);
        v = p.z - floorf(p.z);
    } else if (hit.normal.y) {
        u = p.x - floorf(p.x);
        v = p.z - floorf(p.z);
    } else {
        u = p.x - floorf(p.x);
        v = p.y - floorf(p.y);
    }

    u *= 128, v *= 128;

    int iu = (int)(u * texture->width) % texture->width;
    int iv = (int)(v * texture->height) % texture->height;
    uint texel = texture->pixels[iu + iv * texture->width];
    hit.albedo = RGB8_to_RGBF32(texel);
#endif

#ifdef DEV
    /* Handle special display modes, for debugging */
    if (dev::display_mode != dev::DM::FINAL) {
        if (dev::display_mode != dev::DM::PRIMARY_STEPS) {
            if (hit.depth >= BIG_F32) {
                result.albedo = float4(0.06f); /* 0xFF101010 */
                result.irradiance = -1;
                return result;
            }
        }

        switch (dev::display_mode) {
            case dev::DM::ALBEDO:
                result.albedo = hit.albedo;
                result.irradiance = -1;
                return result;
            case dev::DM::NORMALS:
                result.albedo = float4(fmaxf(hit.normal, 0.0f), 1.0f);
                result.irradiance = -1;
                return result;
            case dev::DM::DEPTH:
                result.albedo =
                    float4(hit.depth / 16.0f, hit.depth / 16.0f, hit.depth / 16.0f, 1.0f);
                result.irradiance = -1;
                return result;
            case dev::DM::PRIMARY_STEPS:
                result.albedo =
                    float4(hit.steps / 64.0f, hit.steps / 64.0f, hit.steps / 64.0f, 1.0f);
                result.irradiance = -1;
                return result;
        }
    }
#endif

    /* Skybox color if the ray missed */
    if (hit.depth >= BIG_F32) {
        result.albedo = skydome.sample_dir(ray.dir);
        result.irradiance = 0;
        return result;
    }

    result.albedo = hit.albedo;

    /* Ambient light */
#if 1

#ifdef DEV
    u32 al_steps = 0;
#endif
    float3 ambient_c = float3(0.0f);
    constexpr f32 SAMPLES = 1;
    for (u32 i = 0; i < SAMPLES; i++) {
        /* Blue noise + R2 (cosine weighted distribution) */
        const float2 raw_noise = bnoise->sample_2d(x, y);
        const f32 quasi_x = fmod(raw_noise.x + R2X_2D * (f64)(frame + i), 1.0f);
        const f32 quasi_y = fmod(raw_noise.y + R2Y_2D * (f64)(frame + i), 1.0f);
        const float3 ambient_dir = cosineweighteddiffusereflection(hit.normal, quasi_x, quasi_y);

        /* Shoot the ambient ray */
        const Ray ambient_ray = Ray(hit_pos, ambient_dir * 16.0f);
#ifdef DEV
        const bool in_shadow = volume->is_occluded(ambient_ray, 1.0f, &al_steps);
#else
        const bool in_shadow = volume->is_occluded(ambient_ray);
#endif

        if (not in_shadow) {
            /* Adjust the samples based on their probability distribution function (PDF) */
            const f32 pdf = dot(ambient_dir, hit.normal) * INVPI; /* (cos(a) / PI) */
            const float3 sample = skydome.sample_dir(ambient_dir);
            // const float3 sample = skydome.sample_voxel_normal(hit.normal);
            ambient_c += clamp_color(sample / pdf, 8.0f);
        }
    }
    /* Divide by the number of samples */
    result.irradiance += (ambient_c * (1.0f / SAMPLES)) * 0.5f;

#ifdef DEV
    /* Handle special display modes, for debugging */
    if (dev::display_mode == dev::DM::AMBIENT_STEPS) {
        result.albedo = float4(al_steps / 64.0f, al_steps / 64.0f, al_steps / 64.0f, 1.0f);
        result.irradiance = -1;
        return result;
    }
#endif

#endif

    /* Directional light */
#ifdef DEV
    u32 dl_steps = 0;
#endif
    for (u32 i = 0; i < 1; i++) {
        /* Jitter light position (soft shadows) */
        const f32 JITTER_DIAMETER = 6.0f / 16.0f;
        const f32 JITTER_RADIUS = JITTER_DIAMETER * 0.5f;

        float3 raw_noise = bnoise->sample_3d(x, y);
        float3 quasi_noise;
        quasi_noise.x = fmod(raw_noise.x + R2X * (f32)frame, 1.0f);
        quasi_noise.y = fmod(raw_noise.y + R2Y * (f32)frame, 1.0f);
        quasi_noise.z = fmod(raw_noise.z + R2Z * (f32)frame, 1.0f);
        const float3 jitter = (quasi_noise * JITTER_DIAMETER - JITTER_RADIUS);
        const float3 sun_dirj = normalize(sun_dir + jitter);

        /* Do nothing if the normal faces away from the light */
        f32 incidence = dot(hit.normal, sun_dirj);
        if (incidence <= 0.0f) continue;

        /* Shoot shadow ray */
        const Ray shadow_ray = Ray(hit_pos, sun_dirj * 16.0f);
#ifdef DEV
        const bool in_shadow = volume->is_occluded(shadow_ray, BIG_F32, &dl_steps);
#else
        const bool in_shadow = volume->is_occluded(shadow_ray);
#endif

        if (not in_shadow) {
            const float3 sun_light = float3(0.5f);
            result.irradiance += sun_light * incidence;
        }
    }

#ifdef DEV
    /* Handle special display modes, for debugging */
    if (dev::display_mode == dev::DM::SECONDARY_STEPS) {
        result.albedo = float4(dl_steps / 64.0f, dl_steps / 64.0f, dl_steps / 64.0f, 1.0f);
        result.irradiance = -1;
        return result;
    }
#endif

    u32 noise_n = 0;
    for (const SphereLight& light : area_lights) {
        /* Compute some quasi random noise (for stochastic sampling) */
        const float3 raw_noise = bnoise->sample_3d(x, y);
        const float3 quasi_noise = make_float3(fmod(raw_noise.x + R2X * (f32)frame, 1.0f),
                                               fmod(raw_noise.y + R2Y * (f32)frame, 1.0f),
                                               fmod(raw_noise.z + R2Z * (f32)frame, 1.0f));

        /* Get the contribution for this light */
        result.irradiance += light.contribution(ray, hit, hit_pos, volume, quasi_noise);
        noise_n++;
    }

#if 0
    for (const LightSource& light : lights) {
#if 1
        /* Jitter light position (soft shadows) */
        const f32 JITTER_DIAMETER = 3.0f / 16.0f;
        const f32 JITTER_RADIUS = JITTER_DIAMETER * 0.5f;

        float3 raw_noise = bnoise.sample_3d(x, y);
        float3 quasi_noise;
        quasi_noise.x = fmod(raw_noise.x + R2X * (f32)frame, 1.0f);
        quasi_noise.y = fmod(raw_noise.y + R2Y * (f32)frame, 1.0f);
        quasi_noise.z = fmod(raw_noise.z + R2Z * (f32)frame, 1.0f);
        const float3 light_pos = light.origin + (quasi_noise * JITTER_DIAMETER - JITTER_RADIUS);
#else
        const float3 light_pos = light.origin;
#endif

        const f32 light_dist = length(light_pos - hit_pos);
        const float3 light_dir = (light_pos - hit_pos) / light_dist;

        /* Do nothing if the normal faces away from the light */
        f32 incidence = dot(hit.normal, light_dir);
        if (incidence <= 0.0f) continue;

        /* Handle aperture of spot light */
        if (light.aperture < 1.0f) {
            const f32 a = (dot(light_dir, light.dir) + 1.0f) * 0.5f;
            if (a > light.aperture) continue;
            incidence *= 1.0f - (a / light.aperture);
        }

        const Ray shadow_ray = Ray(light_pos, hit_pos - light_pos);
        const bool in_shadow = volume->is_occluded(shadow_ray);

        /* Do nothing if the point is in shadow */
        if (in_shadow) continue;
        const f32 sqd = light_dist * light_dist;

        /* Area contribution */
        const float3 area_c = light.light * JITTER_DIAMETER;
        color += hit_color * area_c * incidence / sqd;
    }
#endif

    return result;
}

// TODO: Move this stuff out of "renderer.cpp"
/* Source : <https://www.youtube.com/watch?v=KPoeNZZ6H4s> */
class SecondOrderDyn {
    f32 xp, y, yd;
    f32 k1, k2, k3;

   public:
    SecondOrderDyn(f32 f, const f32 z, const f32 r, const f32& x0) {
        f = fmaxf(f, 0.01f);
        k1 = z / (PI * f);
        k2 = 1 / ((2 * PI * f) * (2 * PI * f));
        k3 = r * z / (2 * PI * f);

        xp = x0;
        y = x0;
        yd = 0;
    }

    f32 update(const f32 dt, const f32 x, const f32 xd) {
        // const f32 xd = (x - xp) / dt;
        const f32 k2_stable = fmaxf(k2, 1.1f * (dt * dt / 4 + dt * k1 / 2));
        y = y + dt * yd;
        yd = yd + dt * (x + k3 * xd - y - k1 * yd) / k2_stable;
        return y;
    }
};

constexpr f32 F = 0.15f;
constexpr f32 Z = 0.5f;
constexpr f32 R = 2.0f;
static SecondOrderDyn dyn[4] = {SecondOrderDyn(F, Z, R, 0), SecondOrderDyn(F, Z, R, 0),
                                SecondOrderDyn(F, Z, R, 0), SecondOrderDyn(F, Z, R, 0)};

/**
 * @brief Called every frame.
 */
void Renderer::tick(f32 dt) {
    frame++;
    if (frame > 120) frame = 0;
    Timer t;

#ifndef PROFILING
#pragma omp parallel for schedule(dynamic)
#endif
    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            Ray ray = camera.get_primary_ray(x, y);

            HitInfo hit; /* Trace the scene */
            const TraceResult r = trace(ray, hit, x, y);

            if (hit.depth >= BIG_F32) {
                const float4 c =
                    aces_approx(r.albedo * insert_accu_raw(x, y, ray, float4(1.0f, hit.depth)));
                albedo_buf[x + y * WIN_WIDTH] = r.albedo * 1.0f;
                screen->pixels[x + y * WIN_WIDTH] = RGBF32_to_RGB8(&c);
                continue;
            }

            if (r.irradiance.x == -1) {
#if DENOISE
                const float4 c = insert_accu_raw(x, y, ray, float4(r.albedo, hit.depth));
#else
                const float4 c = r.albedo;
#endif
                screen->pixels[x + y * WIN_WIDTH] = RGBF32_to_RGB8(&c);
                continue;
            }

            /* Accumulate and reproject */
            const float4 c =
                aces_approx(r.albedo * insert_accu_raw(x, y, ray, float4(r.irradiance, hit.depth)));
            albedo_buf[x + y * WIN_WIDTH] = r.albedo;
            screen->pixels[x + y * WIN_WIDTH] = RGBF32_to_RGB8(&c);
        }
    }

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

    /* Update the camera */
    depth_delta = camera.update(dt);
}

void Renderer::gui(f32 dt) {
    /* Development GUI */
    devgui_stats(dt);
    devgui_control();

#ifdef DEV
    /* Fast mode switch */
    static bool f_down = false;
    if (IsKeyDown(GLFW_KEY_F) && f_down == false) {
        if (dev::display_mode == dev::DM::FINAL) {
            dev::display_mode = dev::DM::ALBEDO;
        } else {
            dev::display_mode = dev::DM::FINAL;
        }
        f_down = true;
    }
    if (!IsKeyDown(GLFW_KEY_F) && f_down == true) {
        f_down = false;
    }

    /* Dev gui switch */
    static bool tilde_down = false;
    if (IsKeyDown(GLFW_KEY_GRAVE_ACCENT) && !tilde_down) {
        dev::hide_devgui = !dev::hide_devgui;
        tilde_down = true;
    }
    if (!IsKeyDown(GLFW_KEY_GRAVE_ACCENT) && tilde_down) {
        tilde_down = false;
    }
#endif

#ifndef PROFILING
    /* TEMP: Draw a cube */
    {
        const float3 point = 0, size = 1;
        static f32 time = 0;
        time += dt;
        const quat rot = quat::from_axis_angle(normalize({1, 1, 0}), time);
        db::draw_obb(point, size, rot);

        /* Cube planes */
        const float3 planes[6] = {
            {0, 1, 0},  /* top */
            {0, -1, 0}, /* bottom */
            {1, 0, 0},  /* right */
            {-1, 0, 0}, /* left */
            {0, 0, 1},  /* forward */
            {0, 0, -1}, /* backward */
        };
        for (u32 i = 0; i < 6; i++) {
            const float3 plane = rot.rotateVector(planes[i] * (size * 0.5f));
            const float3 normal = normalize(plane);
            db::draw_normal(plane, normal, 0xFFFFFF00);
        }
    }
#endif
}

void Renderer::shutdown() {
    /* Save the camera state */
    FILE* f = fopen("camera.bin", "wb");
    fwrite(&camera, 1, sizeof(Camera), f);
    fclose(f);

#if USE_BVH
    delete bvh;
#else
    delete volume;
#endif
    delete bnoise;
}

void Renderer::MouseDown(int button) {
#ifdef DEV
    /* Don't listen to mouse if it's over ImGui windows */
    if (ImGui::GetIO().WantCaptureMouse) return;
#endif

    if (button == 0) {
        // volume->place_voxel(camera.get_primary_ray(mousePos.x, mousePos.y));
    }
}

/**
 * @brief Reproject onto the current frame and accumulate. (with tonemapping)
 * @return The color to display on screen for this pixel.
 */
inline float4 Renderer::insert_accu(const u32 x, const u32 y, const Ray& ray,
                                    const float4& c) const {
    return aces_approx(insert_accu_raw(x, y, ray, c));
}

/**
 * @brief Reproject onto the current frame and accumulate. (without tonemapping)
 * @return The color to display on screen for this pixel.
 */
inline float4 Renderer::insert_accu_raw(const u32 x, const u32 y, const Ray& ray,
                                        const float4& c) const {
#ifdef DEV
    if (not dev::use_projection) return c;
#endif

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
