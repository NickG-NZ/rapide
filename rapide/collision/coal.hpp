#pragma once
#include "rapide/maximal/forces.hpp"
#include <coal/broadphase/broadphase_callbacks.h>
#include <coal/broadphase/broadphase_dynamic_AABB_tree.h>
#include <coal/collision.h>
#include <coal/shape/geometric_shapes.h>
#include <memory>
#include <numeric>
#include <set>

namespace rapide::collision {

using maximal::BodyId;
using maximal::Pose;
using maximal::Vec3;
using ColliderId = int;

/** Normal points A -> B; gap is negative in penetration. Witness points are
 *  in authored body coordinates, while midpoint and normal are world-expressed.
 */
struct Contact {
  ColliderId collider_a, collider_b;
  BodyId body_a, body_b;
  Vec3 local_a, local_b, midpoint, normal;
  double gap;
};

/** Owns Coal objects and a dynamic AABB tree. Model topology is fixed at
 *  construction; shape transforms are refreshed from each queried state.
 *  One world is mutable and must not be queried concurrently. */
class World {
public:
  explicit World(const maximal::Model &model,
                 bool exclude_joint_neighbors = true)
      : model_(model) {
    if (exclude_joint_neighbors) {
      std::vector<int> group(model.bodies.size());
      std::iota(group.begin(), group.end(), 0);
      for (const auto &j : model.joints)
        if (j.type == maximal::JointType::Fixed) {
          const int from = group[j.b], to = group[j.a];
          for (auto &value : group)
            if (value == from)
              value = to;
        }
      std::set<std::pair<int, int>> adjacent;
      for (const auto &j : model.joints)
        adjacent.insert(ordered(group[j.a], group[j.b]));
      for (int a = 0; a < static_cast<int>(group.size()); ++a)
        for (int b = a + 1; b < static_cast<int>(group.size()); ++b)
          if (group[a] == group[b] ||
              adjacent.contains(ordered(group[a], group[b])))
            exclude(a, b);
    }
  }

  World(const World &) = delete;
  World &operator=(const World &) = delete;

  ColliderId add(BodyId body, std::shared_ptr<coal::CollisionGeometry> geometry,
                 const Pose &local = {}) {
    check_body(body);
    maximal::validate_pose(local);
    if (!geometry)
      throw std::invalid_argument("Collider geometry is null");
    auto collider = std::make_unique<Collider>();
    collider->id = static_cast<int>(colliders_.size());
    collider->body = body;
    collider->local = local;
    collider->object = std::make_unique<coal::CollisionObject>(geometry);
    collider->object->setUserData(collider.get());
    manager_.registerObject(collider->object.get());
    colliders_.push_back(std::move(collider));
    manager_.setup();
    return static_cast<int>(colliders_.size() - 1);
  }

  void exclude(BodyId a, BodyId b) {
    check_body(a);
    check_body(b);
    excluded_.insert(ordered(a, b));
  }
  void include(BodyId a, BodyId b) {
    check_body(a);
    check_body(b);
    excluded_.erase(ordered(a, b));
  }

  std::vector<Contact> contacts(const maximal::State &state) {
    model_.validate_state(state);
    for (const auto &collider : colliders_) {
      const auto &s = state.bodies[collider->body];
      const auto q = s.pose.orientation * collider->local.orientation;
      const auto p = maximal::world_point(model_, state, collider->body,
                                          collider->local.position);
      collider->object->setTransform(
          coal::Transform3s(q.toRotationMatrix(), p));
      collider->object->computeAABB();
    }
    manager_.update();
    Callback callback(*this, state);
    manager_.collide(&callback);
    std::sort(callback.result.begin(), callback.result.end(),
              [](const Contact &a, const Contact &b) {
                return std::pair{a.collider_a, a.collider_b} <
                       std::pair{b.collider_a, b.collider_b};
              });
    return std::move(callback.result);
  }

private:
  struct Collider {
    ColliderId id;
    BodyId body;
    Pose local;
    std::unique_ptr<coal::CollisionObject> object;
  };

  struct Callback : coal::CollisionCallBackBase {
    const World &world;
    const maximal::State &state;
    std::vector<Contact> result;
    Callback(const World &w, const maximal::State &s) : world(w), state(s) {}

    bool collide(coal::CollisionObject *first,
                 coal::CollisionObject *second) override {
      auto *a = static_cast<Collider *>(first->getUserData());
      auto *b = static_cast<Collider *>(second->getUserData());
      if (a->id > b->id)
        std::swap(a, b);
      if (a->body == b->body ||
          world.excluded_.contains(ordered(a->body, b->body)) ||
          (!world.model_.bodies[a->body].dynamic() &&
           !world.model_.bodies[b->body].dynamic()))
        return false;
      coal::CollisionRequest request;
      request.enable_contact = true;
      request.num_max_contacts = 8;
      coal::CollisionResult hits;
      coal::collide(a->object.get(), b->object.get(), request, hits);
      for (size_t i = 0; i < hits.numContacts(); ++i) {
        const auto &c = hits.getContact(i);
        const Vec3 pa = local_point(a->body, c.nearest_points[0]);
        const Vec3 pb = local_point(b->body, c.nearest_points[1]);
        if (!pa.allFinite() || !pb.allFinite() || !c.normal.allFinite() ||
            !std::isfinite(c.penetration_depth))
          throw std::runtime_error("Coal returned a nonfinite contact");
        result.push_back({a->id, b->id, a->body, b->body, pa, pb, c.pos,
                          c.normal, c.penetration_depth});
      }
      return false;
    }

    Vec3 local_point(BodyId body, const Vec3 &p) const {
      const auto &s = state.bodies[body];
      return world.model_.bodies[body].com +
             s.pose.orientation.conjugate() * (p - s.pose.position);
    }
  };

  static std::pair<int, int> ordered(int a, int b) {
    return {std::min(a, b), std::max(a, b)};
  }
  void check_body(BodyId body) const {
    if (body < 0 || body >= static_cast<int>(model_.bodies.size()))
      throw std::invalid_argument("Invalid collider body ID");
  }

  maximal::Model model_;
  std::vector<std::unique_ptr<Collider>> colliders_;
  coal::DynamicAABBTreeCollisionManager manager_;
  std::set<std::pair<int, int>> excluded_;
};

struct Compliance {
  double stiffness = 1e5;
  double damping = 1e3;
  double friction = .8;
  double slip_regularization = .01;
};

/** Penalty response per contact point, with regularized sliding friction.
 *  This permits penetration and does not enforce sticking or hard contact. */
inline void add_compliant_forces(const maximal::Model &model,
                                 const maximal::State &state,
                                 const std::vector<Contact> &contacts,
                                 maximal::Forces &forces,
                                 const Compliance &material = {}) {
  for (double x : {material.stiffness, material.damping, material.friction,
                   material.slip_regularization})
    if (!std::isfinite(x) || x < 0)
      throw std::invalid_argument("Invalid contact material");
  if (material.slip_regularization == 0)
    throw std::invalid_argument("Slip regularization must be positive");
  for (const auto &c : contacts) {
    if (c.gap >= 0)
      continue;
    const auto local_midpoint = [&](BodyId body) -> Vec3 {
      const auto &s = state.bodies.at(body);
      return model.bodies.at(body).com +
             s.pose.orientation.conjugate() * (c.midpoint - s.pose.position);
    };
    const Vec3 pa = local_midpoint(c.body_a), pb = local_midpoint(c.body_b);
    const Vec3 relative = maximal::point_velocity(model, state, c.body_b, pb) -
                          maximal::point_velocity(model, state, c.body_a, pa);
    const double vn = relative.dot(c.normal);
    const double fn =
        std::max(0., -material.stiffness * c.gap - material.damping * vn);
    const Vec3 tangent = relative - vn * c.normal;
    const Vec3 f =
        fn * c.normal -
        material.friction * fn * tangent /
            std::sqrt(tangent.squaredNorm() + material.slip_regularization *
                                                  material.slip_regularization);
    maximal::add_force_at_point(model, state, forces, c.body_a, pa, -f);
    maximal::add_force_at_point(model, state, forces, c.body_b, pb, f);
  }
}

} // namespace rapide::collision
