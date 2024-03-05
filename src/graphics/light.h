#pragma once

struct LightSource {
    /* The amount of light expressed with 3 color channels */
    float3 light = float3(1.0f, 1.0f, 1.0f);
    float3 origin, dir;
    f32 aperture = 1.0f; /* 1 means point light, anything below is spotlight */

    LightSource(float3 origin, float3 light) : light(light), origin(origin) {}
    LightSource(float3 origin, float3 dir, f32 aperture, float3 light)
        : light(light), origin(origin), dir(dir), aperture(aperture) {}
};
