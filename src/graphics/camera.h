#pragma once

#include "rays/ray.h"
#include "rays/packet.h"

constexpr float3 UP = float3(0, 1, 0);
constexpr f32 ASPECT_RATIO = (f32)WIN_WIDTH / (f32)WIN_HEIGHT;

struct Camera {
    float3 pos, target;
    float3 tl, tr, bl;

    Camera() {
        pos = float3(0, 0, -2);
        target = float3(0, 0, -1);
        tl = float3(-ASPECT_RATIO, 1, 0);
        tr = float3(ASPECT_RATIO, 1, 0);
        bl = float3(-ASPECT_RATIO, -1, 0);
    }

    /* Get a new primary ray from an X and Y pixel coordinate. */
    inline Ray get_primary_ray(const f32 x, const f32 y) const {
        /* UV coordinates */
        const f32 u = x * (1.0f / WIN_WIDTH), v = y * (1.0f / WIN_HEIGHT);
        const float3 ray_end = tl + u * (tr - tl) + v * (bl - tl);
        return Ray(pos, normalize(ray_end - pos));
    }

    /* Bundle of 4 x 3D vectors. */
    union QuadBundle {
        struct {
            float3 v1;
            float3 v2;
            float3 v3;
            float3 v4;
        };
        struct {
            f128 m1;
            f128 m2;
            f128 m3;
        };
    };

    /* Fast SIMD normalize for 4 x 3D vectors. */
    /* Source: <http://virtuallyrandom.com/part-2-vector3-batch-normalization-fpu-vs-simd/> */
    inline static void normalize_bundle(QuadBundle& bundle) {
        const f128 square_0 = _mm_mul_ps(bundle.m1, bundle.m1);
        const f128 square_1 = _mm_mul_ps(bundle.m2, bundle.m2);
        const f128 square_2 = _mm_mul_ps(bundle.m3, bundle.m3);
        const f128 xpose1_0 = _mm_movelh_ps(square_0, square_1);
        const f128 xpose1_1 = _mm_movehl_ps(square_1, square_0);
        const f128 xpose2_1 = _mm_movelh_ps(xpose1_1, square_2);
        const f128 xpose2_2 = _mm_movehl_ps(square_2, xpose1_1);
        const f128 xpose3_0 = _mm_shuffle_ps(xpose1_0, xpose2_2, _MM_SHUFFLE(2, 0, 2, 0));
        const f128 xpose3_2 = _mm_shuffle_ps(xpose1_0, xpose2_2, _MM_SHUFFLE(3, 1, 3, 1));
        const f128 sum1 = _mm_add_ps(xpose3_0, xpose2_1);
        const f128 lensq = _mm_add_ps(xpose3_2, sum1);
        const f128 rcpsqrt = _mm_rsqrt_ps(lensq);
        const f128 scale_0 = _mm_shuffle_ps(rcpsqrt, rcpsqrt, _MM_SHUFFLE(1, 0, 0, 0));
        const f128 scale_1 = _mm_shuffle_ps(rcpsqrt, rcpsqrt, _MM_SHUFFLE(2, 2, 1, 1));
        const f128 scale_2 = _mm_shuffle_ps(rcpsqrt, rcpsqrt, _MM_SHUFFLE(3, 3, 3, 2));
        bundle.m1 = _mm_mul_ps(bundle.m1, scale_0);
        bundle.m2 = _mm_mul_ps(bundle.m2, scale_1);
        bundle.m3 = _mm_mul_ps(bundle.m3, scale_2);
    }

    /* Get a new primary ray packet from a top left X and Y pixel coordinate. */
    inline RayPacket get_primary_packet(const f32 x, const f32 y) const {
        /* UV coordinates */
        // const f128 xm = _mm_set_ps(x, x + 1, x, x + 1), ym = _mm_set_ps(y, y, y + 1, y + 1);
        // const f128 um = _mm_mul_ps(xm, _mm_set_ps1(1.0f / WIN_WIDTH));
        // const f128 vm = _mm_mul_ps(ym, _mm_set_ps1(1.0f / WIN_HEIGHT));

        /* Bundle of 4 x 3D vectors */
        QuadBundle bundle = {};
        bundle.v1 = get_primary_ray(x, y).dir;
        bundle.v2 = get_primary_ray(x + 1, y).dir;
        bundle.v3 = get_primary_ray(x, y + 1).dir;
        bundle.v4 = get_primary_ray(x + 1, y + 1).dir;
        // bundle.v1 = (tl + um.m128_f32[0] * (tr - tl) + vm.m128_f32[0] * (bl - tl)) - pos;
        // bundle.v2 = (tl + um.m128_f32[1] * (tr - tl) + vm.m128_f32[1] * (bl - tl)) - pos;
        // bundle.v3 = (tl + um.m128_f32[2] * (tr - tl) + vm.m128_f32[2] * (bl - tl)) - pos;
        // bundle.v4 = (tl + um.m128_f32[3] * (tr - tl) + vm.m128_f32[3] * (bl - tl)) - pos;
        //  normalize_bundle(bundle); /* Fast SIMD normalize */

        f128 rd[3] = {}; /* Convert from horizontal to vertical layout */
        rd[0] = _mm_set_ps(bundle.v4.x, bundle.v3.x, bundle.v2.x, bundle.v1.x);
        rd[1] = _mm_set_ps(bundle.v4.y, bundle.v3.y, bundle.v2.y, bundle.v1.y);
        rd[2] = _mm_set_ps(bundle.v4.z, bundle.v3.z, bundle.v2.z, bundle.v1.z);

        /* Shared origin of all four rays */
        const f128 ro[3] = {_mm_set_ps1(pos.x), _mm_set_ps1(pos.y), _mm_set_ps1(pos.z)};

        return RayPacket(ro, rd);
    }

    /* Update the camera and handle inputs. */
    bool update(const f32 t);

   private:
    f32 focal_point = 0.01f;
};
