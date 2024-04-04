#include "materials.h"

#include "graphics/noise/sampler.h"

#include "dev/debug.h"

void eval_glass(MatEval& eval, const Ray& ray, const HitInfo& hit, const Scene& scene,
                const NoiseSampler& noise);

/**
 * @brief Evaluate a material in the scene.
 */
void eval_material(MatEval& eval, const Ray& ray, const HitInfo& hit, const Scene& scene,
                   const NoiseSampler& noise) {
    if (eval.bounces > 8) return;

    /* Intersection point */
    const float3 i = ray.intersection(hit);

    /* Get the hit material row */
    const MaterialRow material = (MaterialRow)floor((hit.material - 1) / 8.0f);

    switch (material) {
        case MaterialRow::GLASS:
            eval_glass(eval, ray, hit, scene, noise);
            break;
        default: /* Diffuse */
            eval.albedo = hit.albedo;
            eval.irradiance = diffuse_light(i, hit.normal, scene, noise);
            if (ray.debug) db::draw_normal(i, hit.normal, 0xFF0000FF);
            break;
    }
}

f32 fresnel_reflect(const f32 n1, const f32 n2, const float3& n, const float3& i, const f32 f0,
                    const f32 f90) {
    /* Schlick aproximation */
    f32 r0 = (n1 - n2) / (n1 + n2);
    r0 *= r0;
    f32 cosX = -dot(n, i);
    if (n1 > n2) {
        const f32 n = n1 / n2;
        const f32 sinT2 = n * n * (1.0 - cosX * cosX);

        /* Total internal reflection */
        if (sinT2 > 1.0) return f90;
        cosX = sqrt(1.0 - sinT2);
    }
    const f32 x = 1.0 - cosX;
    const f32 ret = r0 + (1.0 - r0) * x * x * x * x * x;

    /* Lerp from f0 to f90 */
    return f0 + f90 * ret;
}

/**
 * @brief Evaluate a glass material.
 */
void eval_glass(MatEval& eval, const Ray& ray, const HitInfo& hit, const Scene& scene,
                const NoiseSampler& noise) {
    eval.bounces++;
    constexpr f32 REFRACT_IDX = 1.5f;

    /* Entry ray */
    const float3 entry_dir = refract(hit.normal, ray.dir, 1.0f / REFRACT_IDX);
    Ray i_ray = Ray(ray.intersection(hit), entry_dir);
    i_ray.medium_id = hit.material;

    constexpr u32 MAX_REFLECTIONS = 8;
    const float3 absorption = -(1.0f - float3(hit.albedo));

    f32 mul = 1, absorb_t = 0;
    for (u32 i = 0; i < MAX_REFLECTIONS; ++i) {
        /* Find the next internal exit point */
        const HitInfo i_hit = scene.intersect(i_ray);

        /* Move the ray to the exit point */
        i_ray.origin += i_ray.dir * i_hit.depth;

        /* Beer's law (absorption) */
        absorb_t += i_hit.depth;
        const float3 absorb = exp(absorption * 2.0f * absorb_t);

        /* Reflect / Refract ratio */
        const f32 reflect_mul = fresnel_reflect_prob(REFRACT_IDX, 1.0f, i_ray.dir, i_hit.normal);
        const f32 refract_mul = 1.0f - reflect_mul;

        /* Don't refract if chance is small */
        if (refract_mul < 0.2f) {
            /* Reflect internally, and nudge forward a bit */
            const float3 int_reflect = reflect(i_ray.dir, i_hit.normal);
            i_ray = Ray(i_ray.origin + int_reflect * 0.001f, int_reflect);
            i_ray.medium_id = hit.material;
            continue;
        }

        /* Refraction exit ray (scan ray) */
        const float3 scan_dir = refract(i_hit.normal, i_ray.dir, REFRACT_IDX / 1.0f);
        Ray scan_ray = Ray(i_ray.origin + i_hit.normal * 0.0001f, scan_dir);
        scan_ray.ignore_medium = hit.material;
        const HitInfo scan_info = scene.intersect(scan_ray);

        MatEval scan_eval = eval;
        if (scan_info.no_hit()) {
            /* Skydome was hit */
            scan_eval.albedo = scan_info.albedo;
            scan_eval.irradiance = 1;
        } else {
            /* Evaluate the scanned material */
            eval_material(scan_eval, scan_ray, scan_info, scene, noise);
        }
        /* Apply scan findings */
        eval.albedo += scan_eval.albedo * absorb * refract_mul * mul;
        eval.irradiance += scan_eval.irradiance * absorb * refract_mul * mul;

        /* Don't reflect if chance is small */
        if (reflect_mul < 0.2f || mul < 0.1f) {
            break;
        }

        /* Reflect internally, and nudge forward a bit */
        const float3 int_reflect = reflect(i_ray.dir, i_hit.normal);
        i_ray = Ray(i_ray.origin + int_reflect * 0.001f, int_reflect);
        i_ray.medium_id = hit.material;

        /* Carry the reflection multiplier */
        mul *= reflect_mul;
    }
}

/**
 * @brief Evaluate the irradiance for a given point and normal in the scene.
 */
float3 diffuse_light(const float3& p, const float3& n, const Scene& scene,
                     const NoiseSampler& noise) {
    float3 irradiance = 0;

    /* Evaluate the area lights */
    for (const SphereLight& light : scene.lights) {
        irradiance += light.contribution(p, n, scene, noise);
    }

    /* Pick one at random 50/50 */
    if (RandomFloat() <= 0.5f) {
        /* Evaluate the sun light */
        irradiance += sun_light(p, n, scene, noise) * 2.0f;
    } else {
        /* Evaluate the ambient light */
        irradiance += ambient_light(p, n, scene, noise) * 2.0f;
    }

    return irradiance;
}

/**
 * @brief Evaluate sun light contribution for a given point and normal.
 */
float3 sun_light(const float3& p, const float3& n, const Scene& scene, const NoiseSampler& noise) {
    /* Jitter the sun direction (for soft shadows) */
    const f32 JITTER_INTENSITY = 6.0f / 16.0f, JITTER_HALF = JITTER_INTENSITY * 0.5f;
    const float3 sample = noise.sample_3d() * JITTER_INTENSITY - JITTER_HALF;
    const float3 sun_dir = normalize(scene.sun_dir + sample);

    /* Do nothing if the normal faces away from the light */
    const f32 incidence = dot(n, sun_dir);
    if (incidence <= 0.0f) return 0;

    /* Shoot the shadow ray */
    const Ray shadow_ray = Ray(p, sun_dir);
    const bool in_shadow = scene.is_occluded(shadow_ray, BIG_F32);

    if (not in_shadow) {
        return scene.sun_light * incidence;
    }
    return 0;
}

/**
 * @brief Evaluate ambient light contribution for a given point and normal.
 */
float3 ambient_light(const float3& p, const float3& n, const Scene& scene,
                     const NoiseSampler& noise) {
    /* Maximum check distance for occlusion */
    constexpr f32 MAX_DIST = 1.0f;

    /* Get a cosine weighted diffuse reflection vector */
    const float2 sample = noise.sample_2d();
    const float3 ambient_dir = cos_diffuse_reflect(n, sample.x, sample.y);

    /* Shoot the ambient light ray */
    const Ray ambient_ray = Ray(p, ambient_dir);
    const bool in_shadow = scene.is_occluded(ambient_ray, MAX_DIST);

    if (not in_shadow) {
        /* Adjust the samples based on their probability distribution function (PDF) */
        const f32 pdf = dot(ambient_dir, n) * INVPI; /* (cos(a) / PI) */
        const float3 sample = scene.sample_sky(ambient_ray) * 0.25f;
        return clamp_color(sample / pdf, 8.0f);
    }
    return 0;
}

float fresnel_reflect_prob(const f32 n1, const f32 n2, const float3& n, const float3& incident) {
    // Schlick aproximation
    float r0 = (n1 - n2) / (n1 + n2);
    r0 *= r0;
    float cosX = -dot(n, incident);
    if (n1 > n2) {
        float n = n1 / n2;
        float sinT2 = n * n * (1.0f - cosX * cosX);
        // Total internal reflection
        if (sinT2 > 1.0f) return 1.0f;
        cosX = sqrt(1.0f - sinT2);
    }
    float x = 1.0f - cosX;
    float ret = r0 + (1.0f - r0) * x * x * x * x * x;

    // adjust reflect multiplier for object reflectivity
    ret = (MIN_REFLECT + (1.0f - MIN_REFLECT) * ret);
    return ret;
}

float3 refract(const float3& n, const float3& incident, const f32 eta) {
    const f32 d = dot(n, incident);
    const f32 k = 1.0 - eta * eta * (1.0 - d * d);
    if (k < 0.0)
        return 0;
    else
        return normalize(eta * incident - (eta * d + sqrt(k)) * n);
}
