#pragma once
#include "rapide/aba/tree_data.hpp"
#include "rapide/integration/joint_configuration.hpp"

namespace rapide::aba {

/** Map the model's packed coordinates to the shared per-joint update. */
template <typename S>
inline Eigen::VectorX<S>
integrate_joint_configuration(const TreeData<S> &model,
                              const Eigen::VectorX<S> &configuration,
                              const Eigen::VectorX<S> &velocity, S dt) {
  Eigen::VectorX<S> next = configuration;
  for (int i = 0; i < model.num_joints(); ++i) {
    next.segment(model.q_idx[i], joint_nq(model.joints[i])) =
        integration::integrate_joint_configuration(model.joints[i],
                                                   model.q_i(configuration, i),
                                                   model.v_i(velocity, i), dt);
  }
  return next;
}

} // namespace rapide::aba
