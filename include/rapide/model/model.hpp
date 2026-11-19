#pragma once
#include "rapide/math/se3.hpp"
#include "rapide/model/body.hpp"
#include "rapide/model/joint.hpp"
#include <cassert>
#include <string>
#include <vector>

namespace rapide {

// State<S> — runtime joint configuration and velocity.
// q: concatenated joint configurations (nq scalars total)
// v: concatenated joint velocities    (nv scalars total)
template <typename S> struct State {
  Eigen::VectorX<S> q;
  Eigen::VectorX<S> v;
};

// Model<S> — kinematic tree topology and inertia data.
//
// Body/joint indexing:
//   bodies[0]  = world (fixed, infinite mass — inertia is ignored)
//   bodies[i]  = link i  (i >= 1)
//   joints[i]  = joint connecting bodies[parent[i]] to bodies[i+1]
//   X_tree[i]  = fixed SE3 from parent-body frame to joint predecessor frame
//
// Configuration storage:
//   q_idx[i], v_idx[i] = start offset of joint i's slice in q / v.
template <typename S> struct Model {
  std::vector<Body<S>> bodies;  // bodies[0] = world
  std::vector<Joint<S>> joints; // joints[i] → bodies[i+1]
  std::vector<int> parent;      // parent body index for each joint
  std::vector<SE3<S>>
      X_tree; // fixed transform from parent frame to joint predecessor frame
  std::vector<int> q_idx; // q offset for joint i
  std::vector<int> v_idx; // v offset for joint i
  int nq = 0;
  int nv = 0;

  int num_joints() const { return static_cast<int>(joints.size()); }
  int num_bodies() const { return static_cast<int>(bodies.size()); }

  // Slice of q (generalized co-ords) for joint i.
  Eigen::Ref<const Eigen::VectorX<S>> q_i(const Eigen::VectorX<S> &q,
                                          int i) const {
    return q.segment(q_idx[i], joint_nq(joints[i]));
  }

  // Slice of v (generalzied velocities) for joint i.
  Eigen::Ref<const Eigen::VectorX<S>> v_i(const Eigen::VectorX<S> &v,
                                          int i) const {
    return v.segment(v_idx[i], joint_nv(joints[i]));
  }

  // Build a default-zero state compatible with this model.
  State<S> zero_state() const {
    State<S> s;
    s.q = Eigen::VectorX<S>::Zero(nq);
    s.v = Eigen::VectorX<S>::Zero(nv);
    // Set quaternion w=1 for any FreeJoint so the orientation is identity.
    for (int i = 0; i < num_joints(); ++i) {
      if (std::holds_alternative<FreeJoint<S>>(joints[i])) {
        s.q(q_idx[i] + 3) = S(1); // qw = 1
      }
    }
    return s;
  }
};

// ─── ModelBuilder
// ───────────────────────────────────────────────────────────── Convenience
// builder; avoids constructing the offsets manually.
template <typename S> class ModelBuilder {
public:
  ModelBuilder() {
    // Add world body — inertia is unused.
    model_.bodies.push_back(Body<S>{"world",
                                    {S(0), Eigen::Matrix<S, 3, 1>::Zero(),
                                     Eigen::Matrix<S, 3, 3>::Zero()}});
  }

  // Add a body with a joint connecting it to parent_body_idx.
  // X_parent: fixed SE3 from parent-body frame to joint predecessor frame.
  // Returns the new body's index.
  int add_body(Body<S> body, Joint<S> joint, int parent_body_idx,
               SE3<S> X_parent = SE3<S>::Identity()) {
    assert(parent_body_idx < num_bodies());

    model_.q_idx.push_back(model_.nq);
    model_.v_idx.push_back(model_.nv);
    model_.nq += joint_nq(joint);
    model_.nv += joint_nv(joint);

    model_.parent.push_back(parent_body_idx);
    model_.joints.push_back(std::move(joint));
    model_.X_tree.push_back(X_parent);
    model_.bodies.push_back(std::move(body));

    return num_bodies() - 1;
  }

  int num_bodies() const { return static_cast<int>(model_.bodies.size()); }

  Model<S> build() { return std::move(model_); }

private:
  Model<S> model_;
};

} // namespace rapide
