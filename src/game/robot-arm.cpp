#include "precomp.h"
#include "robot-arm.h"

RobotArm::RobotArm(RobotBone bones[], u32 len) : arm_len(len) { 
	arm = new RobotBone[len];
    for (u32 i = 0; i < len; i++) {
        arm[i] = bones[i];
	}
}
