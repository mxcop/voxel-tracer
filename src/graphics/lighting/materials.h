#pragma once

class NoiseSampler;

/**
 * @brief Material row table.
 */
enum class MaterialRow : u32 {
    GLASS = 0,    /* ID : 0..8 */
    MIRROR = 1,   /* ID : 8..16 */
    METAL = 2,    /* ID : 16..24 */
    UNUSED_3 = 3, /* ID : 24..32 */
    NOT_LIT = 15
};

struct MatEval {
    float3 albedo = 0, irradiance = 0;
    u32 bounces = 0;
};

/**
 * @brief Evaluate a hit material in the scene.
 */
void eval_material(MatEval& eval, const Ray& ray, const HitInfo& hit,
                     const Scene& scene, const NoiseSampler& noise);

/**
 * @brief Find the next ray along the path through the scene.
 */
bool next_path_ray(Ray& ray, const HitInfo& hit);

/**
 * @brief Evaluate the irradiance for a given point and normal in the scene.
 */
float3 diffuse_light(const float3& p, const float3& n, const Scene& scene, const NoiseSampler& noise);

/**
 * @brief Evaluate sun light contribution for a given point and normal.
 */
float3 sun_light(const float3& p, const float3& n, const Scene& scene, const NoiseSampler& noise);

/**
 * @brief Evaluate ambient light contribution for a given point and normal.
 */
float3 ambient_light(const float3& p, const float3& n, const Scene& scene,
                     const NoiseSampler& noise);

constexpr f32 MIN_REFLECT = 0.01f;

// TODO: move this function out
float fresnel_reflect_prob(const f32 n1, const f32 n2, const float3& n, const float3& incident);

float3 refract(const float3& n, const float3& incident, const f32 eta);
