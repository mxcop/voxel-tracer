#pragma once

struct RobotBone {
    float3 offset = 0;
    float3 axis = 0;
    f32 angle = 0;
    f32 limit = 0;

    RobotBone() = default;
    RobotBone(const float3& offset, const float3& axis, const f32 limit, const f32 angle = 0)
        : offset(offset), axis(axis), angle(angle), limit(limit){};
};

struct RobotArm {
    RobotBone* arm = nullptr;
    u32 arm_len = 0;

    RobotArm() = default;
    RobotArm(RobotBone bones[], u32 len);
};
