#include "blue.h"

#include <stb_image.h>

BlueNoise::BlueNoise() {
    sampler_2d = stbi_loadf("assets/noise/LDR_RG01.png", &w_2d, &h_2d, &n_2d, 4);
    assert(sampler_2d && "Failed to load 2D blue noise texture!");
    sampler_3d = stbi_loadf("assets/noise/LDR_RGB1.png", &w_3d, &h_3d, &n_3d, 4);
    assert(sampler_3d && "Failed to load 3D blue noise texture!");
}

BlueNoise::~BlueNoise() {
    delete[] sampler_2d;
    delete[] sampler_3d;
}
