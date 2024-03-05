#pragma once

/* R2 irrationals */
constexpr f32 R2 = 1.22074408460575947536f;
constexpr f32 R2X = 1.0f / R2;
constexpr f32 R2Y = 1.0f / (R2 * R2);
constexpr f32 R2Z = 1.0f / (R2 * R2 * R2);
constexpr f32 R2_2D = 1.32471795724474602596f;
constexpr f32 R2X_2D = 1.0f / R2_2D;
constexpr f32 R2Y_2D = 1.0f / (R2_2D * R2_2D);

/* Textures obtained from : <http://momentsingraphics.de/BlueNoise.html> */
/* Blue noise sampler */
class BlueNoise {
    f32* sampler_2d, *sampler_3d;
    i32 w_2d, w_3d, h_2d, h_3d, n_2d, n_3d;

   public:
    BlueNoise();
    ~BlueNoise();

    BlueNoise(const BlueNoise&) = delete;
    BlueNoise(BlueNoise&&) = default;
    BlueNoise& operator=(const BlueNoise&) = delete;
    BlueNoise& operator=(BlueNoise&&) = default;

    float2 sample_2d(u32 x, u32 y) const {
        x = x % w_2d, y = y % h_2d;
        const u32 i = (x + y * w_2d) * n_2d;
        return float2(sampler_2d[i + 0], sampler_2d[i + 1]);
    }

    float3 sample_3d(u32 x, u32 y) const {
        x = x % w_3d, y = y % h_3d;
        const u32 i = (x + y * w_3d) * n_3d;
        return float3(sampler_3d[i + 0], sampler_3d[i + 1], sampler_3d[i + 2]);
    }
};
