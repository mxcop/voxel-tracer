#include "sphere-light.h"

SphereLight::SphereLight(const float3 origin, const f32 radius, const float3 color, const f32 power)
    : color(color), power(power), radius(radius), origin(origin), aoe_sqr(aprox_aoe_sqr()) {}

float4 SphereLight::contribution(const Ray& pray, const HitInfo& phit, const float3& surface,
                                 const BrickVolume* scene, const float3& noise) const {
    /* Use Monte Carlo integration over the area of the light */
    const f32 diameter = radius * 2.0f;
    const float3 sample_point = origin + (noise * diameter - radius);

    /* Check if the surface position lies within the area of effect */
    const float3 sample_extend = sample_point - surface;
    const f32 dist_sqr = dot(sample_extend, sample_extend);
    if (dist_sqr > aoe_sqr) return float4(0);

    /* Check if the surface normal faces away from the light */
    const f32 sample_dist = sqrtf(dist_sqr);
    const float3 sample_dir = sample_extend / sample_dist;
    const f32 incidence = dot(phit.normal, sample_dir);
    if (incidence <= 0.0f) return float4(0);

    /* Check if there's a clear line of sight between the surface and sample point */
    const Ray shadow_ray = Ray(sample_point, -sample_extend);
    const bool occluded = scene->is_occluded(shadow_ray);
    if (occluded) return float4(0);

    /* Adjust the samples based on their probability distribution function (PDF) */
    const f32 pdf = FOURPI * diameter; /* (4PI * r2) */

    /* Finally calculate the contribution */
    const f32 i = intensity(dist_sqr);
    const float4 sample = color * i * incidence * pdf;
    return phit.albedo * sample;
}

float4 SphereLight::contribution(const Ray& pray, const HitInfo& phit, const float3& surface,
                                 const VoxelVolume* scene, const float3& noise) const {
    /* Use Monte Carlo integration over the area of the light */
    const f32 diameter = radius * 2.0f;
    const float3 sample_point = origin + (noise * diameter - radius);

    /* Check if the surface position lies within the area of effect */
    const float3 sample_extend = sample_point - surface;
    const f32 dist_sqr = dot(sample_extend, sample_extend);
    if (dist_sqr > aoe_sqr) return float4(0);

    /* Check if the surface normal faces away from the light */
    const f32 sample_dist = sqrtf(dist_sqr);
    const float3 sample_dir = sample_extend / sample_dist;
    const f32 incidence = dot(phit.normal, sample_dir);
    if (incidence <= 0.0f) return float4(0);

    /* Check if there's a clear line of sight between the surface and sample point */
    const Ray shadow_ray = Ray(sample_point, -sample_extend);
    const bool occluded = scene->is_occluded(shadow_ray);
    if (occluded) return float4(0);

    /* Adjust the samples based on their probability distribution function (PDF) */
    const f32 pdf = FOURPI * diameter; /* (4PI * r2) */

    /* Finally calculate the contribution */
    const f32 i = intensity(dist_sqr);
    const float4 sample = color * i * incidence * pdf;
    return phit.albedo * sample;
}
