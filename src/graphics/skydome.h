#pragma once

/* HDR Sky dome */
class SkyDome {
    vector<f32> sampler;
    i32 w, h, n;

    /* Source : <https://gist.github.com/volkansalma/2972237> */
    __forceinline f32 atan2_approx(const f32 y, const f32 x) const {
        // http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
        // Volkan SALMA

        const float ONEQTR_PI = PI / 4.0f;
        const float THRQTR_PI = 3.0f * PI / 4.0f;
        float r, angle;
        float abs_y = fabs(y) + 1e-10f;  // kludge to prevent 0/0 condition
        if (x < 0.0f) {
            r = (x + abs_y) / (abs_y - x);
            angle = THRQTR_PI;
        } else {
            r = (x - abs_y) / (x + abs_y);
            angle = ONEQTR_PI;
        }
        angle += (0.1963f * r * r - 0.9817f) * r;
        if (y < 0.0f)
            return (-angle);  // negate if in quad III or IV
        else
            return (angle);
    }

   public:
    SkyDome() = default;
    SkyDome(const char* file_path);

    float3 sample_dir(float3 dir) const {
        const u32 u = floor(w * atan2_approx(dir.z, dir.x) * INV2PI - 0.5f);
        const u32 v = floor(h * acosf(dir.y) * INVPI - 0.5f);
        const u32 i = min(u + v * w, (u32)(w * h)); /* safety clamp */
        // u32 i = (u + v * w) % (w * h);
        return float3(sampler[i * 3], sampler[i * 3 + 1], sampler[i * 3 + 2]);
    }
};
