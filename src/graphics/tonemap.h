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
