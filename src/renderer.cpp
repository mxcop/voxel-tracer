#include <functional>

#include "graphics/lighting/sample.h"
#include "dev/gui.h"
#include "graphics/rays/frustum.h"
#include "graphics/primitives/basic/sphere.h"
#include <graphics/noise/gaussian.h>

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
    albedo_buf = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (albedo_buf) memset(albedo_buf, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    prev_frame = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    if (prev_frame) memset(prev_frame, 0, WIN_WIDTH * WIN_HEIGHT * sizeof(float4));

#if DENOISE
    blur_in = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
    blur_out = (float4*)MALLOC64(WIN_WIDTH * WIN_HEIGHT * sizeof(float4));
#endif

    bnoise = new BlueNoise();
    skydome = SkyDome("assets/kiara_1_dawn_8k.hdr");

    /* Create a voxel volume */
    // volume = new VoxelVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
    // volume = new BrickVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
#if USE_BVH
    // test_plane_vv = new Sphere(
    //     float3(0.0f, 10.0f, 0.0f),
    //     0.25f);  // new OVoxelVolume(float3(0.0f, 10.0f, 0.0f), "assets/vox/crate-16.vox");
    //  shapes[0] = test_plane_vv;
    //  shapes[1] = new OBB(float3(-0.5f, 2.5f, -0.5f), float3(3), float3(0, 0, 1), 1.0f);

    const f32 VOXEL = 1.0f / 20.0f;

    /* Crates */
    shapes[0] = new OVoxelVolume(float3(-3.0f, 2.5f + VOXEL * 6, -0.5f), "assets/vox/crate-10.vox");
    shapes[1] =
        new OVoxelVolume(float3(-3.0f, 2.5f - VOXEL * 3, -0.5f), "assets/vox/crate-16h.vox");
    shapes[2] =
        new OVoxelVolume(float3(-3.0f + VOXEL * 17, 2.5f, -0.5f), "assets/vox/crate-16.vox");

    /* Robot arm */
    arm_vv[0] = new OVoxelVolume(0, "assets/vox/robot-demo/arm-base-azimuth.vox");
    arm_vv[0]->set_pivot(float3(VOXEL * 4.0f, VOXEL * 4.0f, VOXEL * 4.0f));
    arm_vv[1] = new OVoxelVolume(0, "assets/vox/robot-demo/arm-base-altitude.vox");
    arm_vv[1]->set_pivot(float3(VOXEL * 4.0f, VOXEL * 4.0f, VOXEL * 8.0f));
    arm_vv[2] = new OVoxelVolume(0, "assets/vox/robot-demo/arm-altitude.vox");
    arm_vv[2]->set_pivot(float3(VOXEL * 4.0f, VOXEL * 4.0f, VOXEL * 4.0f));
    arm_vv[3] = new OVoxelVolume(0, "assets/vox/robot-demo/arm-hand.vox");
    arm_vv[3]->set_pivot(float3(VOXEL * 8.0f, VOXEL * 0.0f, VOXEL * 8.0f));

    auto* base = new OVoxelVolume(0, "assets/vox/robot-demo/arm-base.vox");
    base->set_pivot(float3(VOXEL * 8.0f, VOXEL * 4.0f, VOXEL * 8.0f));
    base->set_position(float3(0, 3.0f - VOXEL * 6.0f, 0));
    base->set_rotation(quat::from_axis_angle(float3(0, 1, 0), PI));
    shapes[3] = base;
    shapes[4] = arm_vv[0];
    shapes[5] = arm_vv[1];
    shapes[6] = arm_vv[2];
    shapes[7] = arm_vv[3];

    // shapes[6] = new Sphere(float3(-3.0f + VOXEL * 17, 2.5f, -0.5f), 0.5f);
    // shapes[3] = new AABB(float3(-4, -1, -4), float3(4, 0, 4), float3(1));

    bvh = new Bvh(BVH_SHAPES, shapes);

    /* Physics testing */
    //test_obj = world.add_object(PhyObject(new SphereCollider(0, 1.5f), float3(0, 10, 1.0f), 0.01f));
    //test_plane_obj =
    //    world.add_object(PhyObject(new BoxCollider(float3(-4, -1, -4), float3(4, 0, 4)), 0, 0));

    RobotBone bones[] = {RobotBone(float3(0, VOXEL * 3, 0), float3(0, 1, 0), 0.0f),
                         RobotBone(float3(0, 1, 0), float3(0, 0, 1), PI * 0.6f),
                         RobotBone(float3(0, 1.5f - VOXEL * 3, 0), float3(0, 0, 1), PI * 0.6f),
                         RobotBone(float3(0, 1, 0), float3(0, 1, 0), 0.0f)};
    arm = RobotArm(bones, 4);

#else
    volume = new VoxelVolume(float3(0.0f, 0.0f, 0.0f), int3(128, 128, 128));
#endif

    // texture = new Surface("assets/very-serious-test-image.png");
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
        hit_pos = ray.origin + ray.dir * hit.depth /* + hit.normal * 0.000001f*/;
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

inline float _lerp(const float a, const float b, const float t) { return a + t * (b - a); }

/**
 * @brief Called every frame.
 */
void Renderer::tick(f32 dt) {
    frame++;
    if (frame > 120) frame = 0;
    Timer t;

#if USE_BVH
    //world.step(dt);
    //Transform& transform = test_obj->transform;

    //static f32 time = 0;
    //time += dt;
    //if (time < 5) {
    //    transform.position.y = 5.0f;
    //    test_obj->velocity = 0;
    //}

    //if (transform.position.y < -10.0f) {
    //    transform.position.y = 5.0f;
    //    test_obj->velocity = 0;
    //}
    //test_vv->pos = transform.position;
    //test_plane_vv->pos = test_plane_obj->transform.position;

    static f32 angles[4] = {};
    for (u32 i = 0; i < arm.arm_len; i++) {
        RobotBone& bone = arm.arm[i];

        /* Set new goal */
        if (abs(bone.angle - angles[i]) < 0.1f) {
            if (bone.limit > 0)
                angles[i] = (RandomFloat() * 2 - 1) * bone.limit;
            else
                angles[i] = RandomFloat() * TWOPI;
        }
        /* Move to goal */
        else {
            bone.angle = _lerp(bone.angle, angles[i], 0.1f);
            // dt * 4.0f;
        }
    }

    float3 prev_point = float3(0, 3.0f, 0);
    quat rotation;
    //arm_vv[0]->set_rotation(arm.arm[0].axis * arm.arm[0].angle);
    //arm_vv[0]->set_position(arm.arm[0].offset);
    for (u32 i = 0; i < arm.arm_len; i++) {
        // const RobotBone& prev_bone = arm.arm[i - 1];
        const RobotBone& bone = arm.arm[i];

        const quat bone_rot = quat::from_axis_angle(bone.axis, bone.angle);
        rotation = rotation * bone_rot;
        const float3 next_point = prev_point + rotation.rotateVector(bone.offset);

        arm_vv[i]->set_rotation(rotation);
        arm_vv[i]->set_position(prev_point);

        prev_point = next_point;
    }

    // test_plane_vv->set_position(test_plane_obj->transform.position);
    // test_vv->set_position(transform.position);
    bvh->build(BVH_SHAPES, shapes);
    // reset_accu();
#endif

#pragma omp parallel for schedule(dynamic)
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
#if USE_BVH
    const f32 VOXEL = 1.0f / 20;
    static float3 test_pivot = float3(VOXEL * 4.0f, VOXEL * 3.0f, VOXEL * 4.0f);
    static f32 test_angle = 1.0f;
    if (ImGui::DragFloat3("Pivot", &test_pivot.x, VOXEL) ||
        ImGui::DragFloat("Angle", &test_angle, 0.0174533f)) {
        //arm_vv->set_pivot(test_pivot);
        //arm_vv->set_rotation(float3(0, 0, test_angle));
        bvh->build(BVH_SHAPES, shapes);
    }
#endif

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
