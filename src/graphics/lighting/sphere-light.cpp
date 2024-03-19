#include "sphere-light.h"

#include "graphics/noise/sampler.h"

SphereLight::SphereLight(const float3 origin, const f32 radius, const float3 color, const f32 power)
    : color(color), power(power), radius(radius), origin(origin), aoe_sqr(aprox_aoe_sqr()) {}

float3 SphereLight::contribution(const float3& p, const float3& n, const Scene& scene,
                                 const NoiseSampler& noise) const {
    /* Use Monte Carlo integration over the area of the light */
    const f32 diameter = radius * 2.0f;
    const float3 sample = noise.sample_3d();
    const float3 sample_point = origin + (sample * diameter - radius);

    /* Check if the surface position lies within the area of effect */
    const float3 sample_extend = sample_point - p;
    const f32 dist_sqr = dot(sample_extend, sample_extend);
    if (dist_sqr > aoe_sqr) return 0;

    /* Check if the surface normal faces away from the light */
    const f32 sample_dist = sqrtf(dist_sqr);
    const float3 sample_dir = sample_extend / sample_dist;
    const f32 incidence = dot(n, sample_dir);
    if (incidence <= 0.0f) return 0;

    /* Check if there's a clear line of sight between the surface and sample point */
    const Ray shadow_ray = Ray(sample_point, -sample_dir);
    const bool occluded = scene.is_occluded(shadow_ray, sample_dist - 0.01f);
    if (occluded) return 0;

    /* Adjust the samples based on their probability distribution function (PDF) */
    const f32 pdf = FOURPI * diameter; /* (4PI * r2) */

    /* Finally calculate the contribution */
    const f32 i = intensity(dist_sqr);
    return color * i * incidence * pdf;
}
