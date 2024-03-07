#include <functional>

#include "graphics/lighting/sample.h"
#include "graphics/tonemap.h"
#include "dev/gui.h"
#include "graphics/rays/frustum.h"
#include "graphics/primitives/basic/sphere.h"

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

    /* Create the accumulator */
    accu = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (accu) memset(accu, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));

    bnoise = new BlueNoise();
    skydome = SkyDome("assets/kiara_1_dawn_8k.hdr");

    /* Create a voxel volume */
    // volume = new VoxelVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
    // volume = new BrickVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
#if USE_BVH
    //constexpr u32 SIZE = 16;
    //u32 seed = 47324894723;
    //Traceable** boxes = new Traceable* [SIZE * SIZE + 3] {};
    //for (u32 y = 0; y < SIZE; y++) {
    //    for (u32 x = 0; x < SIZE; x++) {
    //        const f32 xm = x * 8, ym = 0, zm = y * 8;
    //        const u32 i = (y * SIZE) + (0 * SIZE) + x;

    //        const f32 r = RandomFloat(seed);
    //        const u32 c = HSBtoRGB((RandomFloat(seed) * 2 - 1) * 20.0f, 1.0f, 1.0f);
    //        boxes[i] = new AABB(float3(xm, ym, zm),
    //                            float3(xm + 8, ym + r * 8.0f, zm + 8), RGB8_to_RGBF32(c));
    //    }
    //}
    //boxes[SIZE * SIZE] = new Sphere(float3(32 - 4, 9, 40 - 4), 2.5f);
    //boxes[SIZE * SIZE + 1] = new Sphere(float3(16 - 4, 6.5f, 40 - 4), 1.5f);
    //boxes[SIZE * SIZE + 2] = new Sphere(float3(32 - 4, 8, 24 - 4), 2.0f);
    //bvh = new Bvh(SIZE * SIZE + 3, boxes);

    shapes[0] = new AABB(float3(0), float3(1), float3(1));
    shapes[1] = new OBB(float3(-0.5f, 2.5f, -0.5f), float3(3), float3(0, 0, 1), 1.0f);
    test_vv = new OVoxelVolume(float3(2.0f, 2.5f, -0.5f), int3(32), 8);
    test_vv->set_rotation(normalize(RandomFloat3()), RandomFloat() * TWOPI);
    shapes[2] = test_vv;

    bvh = new Bvh(3, shapes);
#else
    volume = new VoxelVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
#endif

    // texture = new Surface("assets/very-serious-test-image.png");
}

u32 Renderer::trace(Ray& ray, const u32 x, const u32 y) const {
#if USE_BVH
    const Bvh* volume = this->bvh;
#endif
    HitInfo hit = volume->intersect(ray);

    float4 color = float4(0);
    float3 hit_pos = ray.origin + ray.dir * hit.depth + hit.normal * 0.00001f;

    /* Reflections */
#if 1
    constexpr u32 MAX_BOUNCES = 8;
    for (u32 i = 0; i < MAX_BOUNCES && hit.albedo.w > 0.0f; ++i) {
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
        hit_pos = ray.origin + ray.dir * hit.depth/* + hit.normal * 0.000001f*/;
    }
#endif

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
            if (hit.depth >= BIG_F32) return 0xFF101010;
        }

        float4 dc;
        switch (dev::display_mode) {
            case dev::DM::ALBEDO:
                return RGBF32_to_RGB8(&hit.albedo);
            case dev::DM::NORMALS:
                // dc = float4(hit.normal + 1.0f * 0.5f, 1.0f);
                dc = float4(fmaxf(hit.normal, 0.0f), 1.0f);
                return RGBF32_to_RGB8(&dc);
            case dev::DM::DEPTH:
                dc = float4(hit.depth / 16.0f, hit.depth / 16.0f, hit.depth / 16.0f, 1.0f);
                return RGBF32_to_RGB8(&dc);
            case dev::DM::PRIMARY_STEPS:
                dc = float4(hit.steps / 64.0f, hit.steps / 64.0f, hit.steps / 64.0f, 1.0f);
                return RGBF32_to_RGB8(&dc);
        }
    }
#endif

    /* Skybox color if the ray missed */
    if (hit.depth >= BIG_F32) {
        color = skydome.sample_dir(ray.dir);

        /* Update accumulator */
        accu[x + y * WIN_WIDTH] += aces_approx(color);
        color = accu[x + y * WIN_WIDTH] / (f32)accu_len;

        return RGBF32_to_RGB8(&color);
    }

    const float3 hit_color = float3(hit.albedo);

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
        const f32 quasi_x = fmod(raw_noise.x + R2X_2D * (frame + i), 1.0f);
        const f32 quasi_y = fmod(raw_noise.y + R2Y_2D * (frame + i), 1.0f);
        // const float3 ambient_dir = sample_hemisphere_weighted(quasi_x, quasi_y, hit.normal);
        const float3 ambient_dir = cosineweighteddiffusereflection(hit.normal, quasi_x, quasi_y);

        /* Shoot the ambient ray */
        const Ray ambient_ray = Ray(hit_pos, ambient_dir * 16.0f);
#ifdef DEV
        const bool in_shadow = volume->is_occluded(ambient_ray, BIG_F32, &al_steps);
#else
        const bool in_shadow = volume->is_occluded(ambient_ray);
#endif

        if (not in_shadow) {
            /* Adjust the samples based on their probability distribution function (PDF) */
            const f32 pdf = dot(ambient_dir, hit.normal) * INVPI; /* (cos(a) / PI) */
            const float3 sample = skydome.sample_dir(ambient_dir);
            ambient_c += (sample / pdf);
        }
    }
    /* Divide by the number of samples */
    color += hit_color * (ambient_c * (1.0f / SAMPLES));

#ifdef DEV
    /* Handle special display modes, for debugging */
    if (dev::display_mode == dev::DM::AMBIENT_STEPS) {
        float4 dc;
        dc = float4(al_steps / 64.0f, al_steps / 64.0f, al_steps / 64.0f, 1.0f);

        /* Update accumulator */
        accu[x + y * WIN_WIDTH] += dc;
        dc = accu[x + y * WIN_WIDTH] / (f32)accu_len;

        return RGBF32_to_RGB8(&dc);
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
            const float3 sun_light = float3(2.5f, 2.5f, 2.5f);
            color += hit_color * sun_light * incidence;
        }
    }

#ifdef DEV
    /* Handle special display modes, for debugging */
    if (dev::display_mode == dev::DM::SECONDARY_STEPS) {
        float4 dc;
        dc = float4(dl_steps / 64.0f, dl_steps / 64.0f, dl_steps / 64.0f, 1.0f);

        /* Update accumulator */
        accu[x + y * WIN_WIDTH] += dc;
        dc = accu[x + y * WIN_WIDTH] / (f32)accu_len;

        return RGBF32_to_RGB8(&dc);
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
        color += light.contribution(ray, hit, hit_pos, volume, quasi_noise);
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

    /* Update accumulator */
    accu[x + y * WIN_WIDTH] += aces_approx(color);
    color = accu[x + y * WIN_WIDTH] / (f32)accu_len;

    return RGBF32_to_RGB8(&color);
}

void Renderer::tick(f32 dt) {
    frame++;
    if (frame > 120) frame = 0;
    Timer t;

#if 0
    // #pragma omp parallel for schedule(dynamic)
    //     for (i32 y = 0; y < WIN_HEIGHT / 12; y++) {
    // #pragma omp parallel for schedule(dynamic)
    //         for (i32 x = 0; x < WIN_WIDTH / 12; x++) {
    //             const u32 xs = x * 12, ys = y * 12;
    //             const Ray ray0 = camera.get_primary_ray(xs, ys);
    //             const Ray ray1 = camera.get_primary_ray(xs + 11, ys);
    //             const Ray ray2 = camera.get_primary_ray(xs, ys + 11);
    //             const Ray ray3 = camera.get_primary_ray(xs + 11, ys + 11);
    //             const Frustum frustum(ray0.origin, ray0.dir, ray1.dir, ray2.dir, ray3.dir,
    //             100.0f);
    //
    //             if (frustum.intersect_unitcube()) {
    //                 for (u32 i = 0; i < 12; i++) {
    //                     for (u32 j = 0; j < 12; j++) {
    //                         const u32 xj = j + xs, yi = i + ys;
    //                         const Ray ray = camera.get_primary_ray(xj, yi);
    //                         const f32 r = bvh->intersect(ray);
    //                         const float4 dc = float4(r / 16.0f, r / 16.0f, r / 16.0f, 1.0f);
    //                         const u32 color = RGBF32_to_RGB8(&dc);
    //                         screen->pixels[xj + yi * WIN_WIDTH] = color;
    //
    //                         // screen->pixels[(j + xs) + (i + ys) * WIN_WIDTH] |= 0xFF00FF00;
    //                     }
    //                 }
    //             }
    //         }
    //     }

#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            Ray ray = camera.get_primary_ray(x, y);
            const HitInfo hit = bvh->intersect(ray);
            // const float4 dc = float4(hit.depth / 64.0f, hit.depth / 64.0f, hit.depth / 64.0f, 1.0f);
            float4 dc;
            switch (dev::display_mode) {
                case dev::DM::ALBEDO:
                    dc = hit.albedo;
                    break;
                case dev::DM::NORMALS:
                    dc = float4(fmaxf(hit.normal, 0.0f), 1.0f);
                    break;
                case dev::DM::DEPTH:
                    dc = float4(hit.depth / 64.0f, hit.depth / 64.0f, hit.depth / 64.0f, 1.0f);
                    break;
                case dev::DM::PRIMARY_STEPS:
                    dc = float4(hit.steps / 64.0f, hit.steps / 64.0f, hit.steps / 64.0f, 1.0f);
                    break;
                default:
                    dc = float4(hit.depth / 64.0f, hit.depth / 64.0f, hit.depth / 64.0f, 1.0f);
                    break;
            }

            u32 color = RGBF32_to_RGB8(&dc);
            screen->pixels[x + y * WIN_WIDTH] = color;
        }
    }
#elif 1
#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; ++y) {
        for (i32 x = 0; x < WIN_WIDTH; ++x) {
            Ray ray = camera.get_primary_ray(x, y);
            u32 color = trace(ray, x, y);
            screen->pixels[x + y * WIN_WIDTH] = color;
        }
    }

    // for (i32 y = 0; y < WIN_HEIGHT / 12; y++) {
    //     for (i32 x = 0; x < WIN_WIDTH / 12; x++) {
    //         u32 xs = x * 12, ys = y * 12;
    //         Ray ray0 = camera.get_primary_ray(xs, ys);
    //         Ray ray1 = camera.get_primary_ray(xs + 11, ys);
    //         Ray ray2 = camera.get_primary_ray(xs, ys + 11);
    //         Ray ray3 = camera.get_primary_ray(xs + 11, ys + 11);
    //         Frustum frustum(ray0.origin, ray0.dir, ray1.dir, ray2.dir, ray3.dir, 100.0f);

    //        if (frustum.intersect_unitcube()) {
    //            for (u32 i = 0; i < 12; i++) {
    //                for (u32 j = 0; j < 12; j++) {
    //                    //u32 xj = j + xs, yi = i + ys;
    //                    //Ray ray = camera.get_primary_ray(xj, yi);
    //                    //u32 color = trace(ray, xj, yi);
    //                    //screen->pixels[xj + yi * WIN_WIDTH] = color;

    //                    screen->pixels[(j + xs) + (i + ys) * WIN_WIDTH] |= 0xFF00FF00;
    //                }
    //            }
    //        }
    //    }
    //}
#else
    // 12.6M rays/s
#pragma omp parallel for schedule(dynamic)
    for (i32 y = 0; y < WIN_HEIGHT; y += 2) {
        for (i32 x = 0; x < WIN_WIDTH; x += 2) {
            const RayPacket128 packet = camera.get_primary_packet(x, y);
            const PacketHitInfo hit = volume->intersect(packet);

            for (u32 v = 0; v < 2; ++v) {
                for (u32 u = 0; u < 2; ++u) {
                    const f32 depth = hit.depth.m128_f32[u + v * 2];
                    /*if (depth == 0.0f) {
                        screen->pixels[(x + u) + (y + v) * WIN_WIDTH] = 0xFFFF1010;
                        continue;
                    }*/
                    if (depth >= BIG_F32) {
                        screen->pixels[(x + u) + (y + v) * WIN_WIDTH] = 0xFF101010;
                        continue;
                    }
                    const u32 cd = fminf(depth / 32.0f, 1.0f) * 0xFF;
                    // const u32 cd = (hit.steps / 256.0f) * 0xFF;
                    const u32 color = (cd << 0) | (cd << 8) | (cd << 16) | 0xFF000000;
                    screen->pixels[(x + u) + (y + v) * WIN_WIDTH] = color;
                }
            }
        }
    }
#endif

    accu_len++;
#ifdef DEV
    dev::frame_time = t.elapsed();
#endif

#if USE_BVH
    //test_vv->set_rotation(normalize(float3(1, 1, 0)), (frame * 3) * 0.0174533f);
    //bvh->build(3, shapes);
#endif

    /* Update the camera */
    if (camera.update(dt)) {
        accu_len = 1u;
        memset(accu, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    }
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
        reset_accu();
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
        reset_accu();
    }
}
