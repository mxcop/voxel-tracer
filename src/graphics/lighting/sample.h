#pragma once

/* Source : <https://github.com/tqjxlm/Monte-Carlo-Ray-Tracer> */
/* Cosine weighted hemisphere sample around normal */
static float3 sample_hemisphere_weighted(const f32 r1, const f32 r2, const float3& n) {
    f32 theta = acos(sqrt(1.0f - r1));
    f32 phi = 2.0f * PI * r2;
    f32 xs = sinf(theta) * cosf(phi);
    f32 ys = cosf(theta);
    f32 zs = sinf(theta) * sinf(phi);
    float3 h = n;

    if ((abs(h.x) <= abs(h.y)) && (abs(h.x) <= abs(h.z))) {
        h.x = 1.0;
    } else if ((abs(h.y) <= abs(h.x)) && (abs(h.y) <= abs(h.z))) {
        h.y = 1.0;
    } else {
        h.z = 1.0;
    }

    float3 x = normalize(cross(h, n));
    float3 z = normalize(cross(x, n));

    return normalize(xs * x + ys * n + zs * z);
}

/* Uniform hemisphere sample around normal */
static float3 sample_hemisphere_uniform(const float3& n) {
    float3 R;
    do {
        R = float3(RandomFloat() * 2 - 1, RandomFloat() * 2 - 1, RandomFloat() * 2 - 1);
    } while (dot(R, R) > 1.0f);
    R = normalize(R);
    if (dot(n, R) < 0.0f) {
        return -R;
    } else {
        return R;
    }
}
