#pragma once

#include "engine/datastruct/pool.h"

#include "object.h"

constexpr f32 GRAVITY = -9.81f;

/**
 * @brief Physics world instance.
 */
class PhyWorld {
    Pool<PhyObject> object_pool;

   public:
    PhyWorld();

    void step(const f32 dt);

    /**
     * @brief Add a new physics object to the world.
     * @return Pointer to the new object. (nullptr if pool was full)
     */
    PhyObject* add_object(const PhyObject obj);

    /**
     * @brief Remove an existing physics object from the world.
     * @return True if it was removed succesfully.
     */
    bool remove_object(const PhyObject* obj);
};
