#pragma once

class BlueNoise;

class NoiseSampler {
    const BlueNoise* sampler = nullptr;
    const u32 x, y, frame;

   public:
    NoiseSampler() = delete;
    /**
     * @brief Create a new noise sampler instance.
     *
     * @param sampler The blue noise sampler to use internally.
     * @param x X axis offset.
     * @param y Y axis offset.
     * @param frame Frame offset. (keep between 0-120 to avoid floating point errors)
     */
    NoiseSampler(const BlueNoise* sampler, const u32 x, const u32 y, const u32 frame);

    /** @brief Sample a 3D noise vector. */
    inline float3 sample_3d(const f32 offset = 0.0f) const {
        const float3 noise = sampler->sample_3d(x, y);
        const f32 frame_offset = (f32)frame + offset;
        return {fmod(noise.x + (f32)R2X * frame_offset, 1.0f),
                fmod(noise.y + (f32)R2Y * frame_offset, 1.0f),
                fmod(noise.z + (f32)R2Z * frame_offset, 1.0f)};
    }

    /** @brief Sample a 2D noise vector. */
    inline float2 sample_2d(const f32 offset = 0.0f) const {
        const float2 noise = sampler->sample_2d(x, y);
        const f32 frame_offset = (f32)frame + offset;
        return {fmod(noise.x + (f32)R2X_2D * frame_offset, 1.0f),
                fmod(noise.y + (f32)R2Y_2D * frame_offset, 1.0f)};
    }
};
