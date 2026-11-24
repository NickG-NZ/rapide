#pragma once
#include "rapide/math/se3.hpp"
#include "rapide/model/model.hpp"

namespace rapide::examples {

/** Attach a serial arm to an authored body; accepted by either backend. */
template <typename S>
void attach_arm(ModelBuilder &builder, int parent_body_idx, int n_dof,
                S link_length = S(.4), S link_mass = S(2),
                const SE3<S> &mount = SE3<S>::Identity()) {
  auto model = builder.build();
  const auto &parent = model.initial_state.bodies.at(parent_body_idx);
  Pose origin{parent.pose.position -
                  parent.pose.orientation * model.bodies[parent_body_idx].com,
              parent.pose.orientation};
  Pose local{mount.p.template cast<double>(),
             Eigen::Quaterniond(mount.R.template cast<double>())};
  int predecessor = parent_body_idx;
  for (int i = 0; i < n_dof; ++i) {
    const Vec3 axis = i % 3 == 1 ? Vec3::UnitY() : Vec3::UnitZ();
    const auto joint_rotation =
        Eigen::Quaterniond::FromTwoVectors(Vec3::UnitZ(), axis);
    const Pose world{origin.position + origin.orientation * local.position,
                     origin.orientation * local.orientation};
    const int body =
        builder.add_body("arm_link_" + std::to_string(i), double(link_mass),
                         Vec3{.01, .01, .001}.asDiagonal(), world,
                         {double(link_length) / 2, 0, 0});
    builder.add_joint(JointType::Revolute, predecessor, body,
                      {local.position, local.orientation * joint_rotation},
                      {{0, 0, 0}, joint_rotation});
    predecessor = body;
    origin = world;
    local = {{double(link_length), 0, 0}};
  }
  const Pose payload_pose{origin.position + origin.orientation * local.position,
                          origin.orientation * local.orientation};
  const auto payload =
      builder.add_body("arm_ee", 5, .001 * Mat3::Identity(), payload_pose);
  const auto wrist =
      Eigen::Quaterniond::FromTwoVectors(Vec3::UnitZ(), Vec3::UnitY());
  builder.add_joint(JointType::Revolute, predecessor, payload,
                    {local.position, local.orientation * wrist},
                    {{0, 0, 0}, wrist});
}

} // namespace rapide::examples
