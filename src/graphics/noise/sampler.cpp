#include "sampler.h"

NoiseSampler::NoiseSampler(const BlueNoise* sampler, const u32 x, const u32 y, const u32 frame)
    : sampler(sampler), x(x), y(y), frame(frame) {}
