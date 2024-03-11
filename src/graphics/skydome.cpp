#include "skydome.h"

#include <stb_image.h>

SkyDome::SkyDome(const char* file_path) {
    f32* data = stbi_loadf(file_path, &w, &h, &n, 0);
    sampler = vector<f32>(data, data + (w * h * n));

    for (f32& sample : sampler) {
        sample = sqrtf(sample);
    }

    /* Integrate voxel normal colors */
    const u32 SAMPLES = 128;
    for (u32 i = 0; i < 3; i++) {
        float3 normal = float3(0);
        normal[i] = 1;

        float3 accu = float3(0);
        for (u32 s = 0; s < SAMPLES; s++) {
            const float3 dir = cosineweighteddiffusereflection(
                normal, RandomFloat() * 0.95f + 0.025f, RandomFloat() * 0.95f + 0.025f);
            accu += sample_dir(dir);
        }
        int_samples[i] = accu / SAMPLES;
    }
    for (u32 i = 0; i < 3; i++) {
        float3 normal = float3(0);
        normal[i] = -1;

        float3 accu = float3(0);
        for (u32 s = 0; s < SAMPLES; s++) {
            const float3 dir = cosineweighteddiffusereflection(
                normal, RandomFloat() * 0.95f + 0.025f, RandomFloat() * 0.95f + 0.025f);
            accu += sample_dir(dir);
        }
        int_samples[i + 3] = accu / SAMPLES;
    }
}
