#include "precomp.h"
#include "world.h"

PhyWorld::PhyWorld() : object_pool(Pool<PhyObject>(64)) {}

void PhyWorld::step(const f32 dt) {
    for (PhyObject& obj : object_pool) {
        if (obj.mass == 0) continue;

        /* Apply gravity */
        obj.force += obj.mass * float3(0, GRAVITY, 0);

        /* Apply force and velocity */
        obj.velocity += obj.force / obj.mass * dt;
        obj.position += obj.velocity * dt;

        obj.force = 0;
    }
}

PhyObject* PhyWorld::add_object(const PhyObject obj) {
    const i32 i = object_pool.add(obj);
    if (i == -1) return nullptr;

    /* Set the objects unique ID */
    PhyObject* object = &object_pool.get(i);
    object->uid = i;
    return object;
}

bool PhyWorld::remove_object(const PhyObject* obj) { return object_pool.remove(obj->uid); }
