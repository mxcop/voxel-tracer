#pragma once

inline float3 reinhard(const float3& v) { return v / (1.0f + v); }

inline float3 reinhard_extended(const float3& v, f32 max_white) {
    const float3 numerator = v * (1.0f + (v / float3(max_white * max_white)));
    return numerator / (1.0f + v);
}

inline float3 aces_approx(float3 v) {
    v *= 0.6f;
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((v * (a * v + b)) / (v * (c * v + d) + e), 0.0f, 1.0f);
}

inline float3 uncharted2(float3 v) {
    const float gamma = 2.4f;
    float A = 0.15f;
    float B = 0.50f;
    float C = 0.10f;
    float D = 0.20f;
    float E = 0.02f;
    float F = 0.30f;
    float W = 11.2f;
    float exposure = 2.0f;
    v *= exposure;
    v = ((v * (A * v + C * B) + D * E) / (v * (A * v + B) + D * F)) - E / F;
    float white = ((W * (A * W + C * B) + D * E) / (W * (A * W + B) + D * F)) - E / F;
    v /= white;
    float3 igamma = float3(1.0f / gamma);
    v.x = powf(v.x, igamma.x);
    v.y = powf(v.y, igamma.y);
    v.z = powf(v.z, igamma.z);
    return v;
}
