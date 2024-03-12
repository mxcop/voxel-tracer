#include "precomp.h"
#include "colliders.h"

float3 SphereCollider::furthest_point(const Transform& t, const float3& dir) const {
    return (t.position + center) + normalize(dir) * radius;
}

float3 BoxCollider::furthest_point(const Transform& t, const float3& dir) const { 
	float3 maxp;
    f32 maxd = -BIG_F32;

    for (u32 i = 0; i < 8; i++) {
        const f32 x = (i & 1) ? max.x : min.x;
        const f32 y = (i & 2) ? max.y : min.y;
        const f32 z = (i & 4) ? max.z : min.z;
        const float3 corner = t.position + float3(x, y, z);

        const f32 d = dot(corner, dir);
        if (d > maxd) {
            maxd = d;
            maxp = corner;
        }
    }

    return maxp;
}
