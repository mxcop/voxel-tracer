#pragma once

/**
 * @brief Spherical light source, radiating light into the scene from a sphere.
 */
class SphereLight {
    /* Color of the light emitted. (RGB) */
    float3 color;
    /* Power the light emits into the scene. */
    f32 power;
    /* Radius of the light sphere. */
    const f32 radius;
    /* Center of the light sphere. */
    float3 origin;
    /* Area of effect of the light. (distance squared) */
    const f32 aoe_sqr;

    /**
     * @brief Calculate the aproximate area of effect of the light. (distance squared)
     *
     * @param power The total power of the light.
     * @return Area of effect of the light as distance squared.
     */
    inline f32 aprox_aoe_sqr() const {
        constexpr f32 THRESHOLD = 1.0f;

        /* Use the inverse-square law to calculate the area of effect distance (squared!) */
        return power * (1.0f / (THRESHOLD * FOURPI));
        /* D = P / (I * 4PI) */
    }

    /**
     * @brief Calculate the intensity of the light at a certain distance. (distance squared)
     *
     * @param dist_sqr The distance squared. (saves a square root)
     * @return Intensity of the light at distance squared.
     */
    inline f32 intensity(const f32 dist_sqr) const {
        /* Use the inverse-square law to calculate the intensity of the light. */
        return power / (FOURPI * dist_sqr);
        /* I = P / (4PI * D) */
    }

   public:
    SphereLight() = delete;
    SphereLight(const float3 origin, const f32 radius, const float3 color, const f32 power);

    float4 contribution(const Ray& pray, const HitInfo& phit, const float3& surface,
                        const BrickVolume* scene, const float3& noise) const;
    float4 contribution(const Ray& pray, const HitInfo& phit, const float3& surface,
                        const VoxelVolume* scene, const float3& noise) const;
};
