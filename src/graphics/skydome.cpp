#include "skydome.h"

#include <stb_image.h>

SkyDome::SkyDome(const char* file_path) {
    f32* data = stbi_loadf(file_path, &w, &h, &n, 0);
    sampler = vector<f32>(data, data + (w * h * n));

    for (f32& sample : sampler) {
        sample = sqrtf(sample);
    }
}
