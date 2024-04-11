#include "world.h"

#include "collision/collision.h"

PhyWorld::PhyWorld() : object_pool(Pool<PhyObject>(64)) {}

void PhyWorld::step(const f32 dt) {
    resolve(dt);

    for (PhyObject& obj : object_pool) {
        if (obj.mass == 0) continue;

        /* Apply gravity */
        obj.force += obj.mass * float3(0, GRAVITY, 0);

        /* Apply force and velocity */
        obj.velocity += obj.force / obj.mass * dt;
        obj.transform.position += obj.velocity * dt;

        obj.force = 0;
    }
}

void PhyWorld::resolve(const f32) {
    /* Collision buffer */
    vector<Collision> c_buffer = {};
    c_buffer.reserve(object_pool.get_size());

    for (PhyObject& a : object_pool) {
        for (PhyObject& b : object_pool) {
            /* Don't collide with yourself */
            if (&a == &b) break;

            /* Check if any colliders are missing */
            if (!a.collider || !b.collider) continue;

            CollisionPoints points =
                collision_test(a.collider, &a.transform, b.collider, &b.transform);
            if (points.collision) c_buffer.emplace_back(&a, &b, points);
        }
    }

    for (Collision& c : c_buffer) {
        /* Testing collision! */
        //if (c.a->mass) {
        //    c.a->transform.position.y = 10.0f;
        //}
        //if (c.b->mass) {
        //    c.b->transform.position.y = 10.0f;
        //}
        PhyObject* a_body = c.a;
        PhyObject* b_body = c.b;

        //a_body->velocity = 0;
        //b_body->velocity = 0;

        //const f32 a_static = (int)(a_body->mass == 0);
        //const f32 b_static = (int)(b_body->mass == 0);

        //const float3 resolution =
        //    c.points.normal * c.points.dist / fmaxf(1, a_static + b_static);

        //a_body->transform.position += resolution * (1.0f - a_static);
        //b_body->transform.position -= resolution * (1.0f - b_static);

        a_body->velocity = 0;
        b_body->velocity = 0;
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
