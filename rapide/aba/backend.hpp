#pragma once
#include "rapide/aba/configuration.hpp"
#include "rapide/aba/dynamics.hpp"
#include "rapide/aba/kinematics.hpp"
#include "rapide/backend.hpp"
#include "rapide/integration/euler.hpp"
#include <numeric>
#include <queue>

namespace rapide::aba {

/** Compiles an acyclic authored graph into ABA frames and joint coordinates. */
class ABABackend final : public SimulationBackend {
public:
  /** Reject cycles (including paths between static mounts), then root each
   *  component at a static mount or an implicit floating joint. Compile bodies
   *  into parent-first order with inertia expressed in child attachment frames,
   *  and extract joint coordinates from the supplied initial body state. */
  void setup(const rapide::Model &model,
             const SolverOptions &options) override {
    options.validate();
    model_ = model;
    tolerance_ =
        std::max({options.position_tolerance, options.orientation_tolerance,
                  options.velocity_tolerance, 1e-8});
    const int count = static_cast<int>(model.bodies.size());
    if (count == 0 || model.bodies[0].dynamic()) {
      throw std::invalid_argument(
          "Model needs a static world body at index zero");
    }
    std::vector<int> component(count);
    std::iota(component.begin(), component.end(), 0);
    const auto find = [&](int body) {
      while (component[body] != body) {
        body = component[body];
      }
      return body;
    };
    for (int i = 0; i < count; ++i) {
      if (!model.bodies[i].dynamic()) {
        component[i] = 0;
      }
    }
    std::vector<std::vector<int>> adjacent(count);
    for (int i = 0; i < static_cast<int>(model.joints.size()); ++i) {
      const auto &joint = model.joints[i];
      joint.validate();
      if (joint.a >= count || joint.b >= count) {
        throw std::invalid_argument("Joint body index outside model");
      }
      const int a = find(joint.a), b = find(joint.b);
      if (a == b) {
        throw std::invalid_argument(
            "ABA does not support kinematic loops; use the maximal backend");
      }
      component[b] = a;
      adjacent[joint.a].push_back(i);
      adjacent[joint.b].push_back(i);
    }
    for (const auto &joint : model.joints) {
      if (joint.type != JointType::Fixed && joint.type != JointType::Revolute &&
          joint.type != JointType::Prismatic) {
        throw std::invalid_argument(
            "ABA supports fixed, revolute and prismatic attachments; floating "
            "roots are implicit");
      }
    }

    TreeBuilder<double> builder;
    body_to_node_.assign(count, 0);
    node_to_body_ = {0};
    frames_.assign(count, SE3<double>::Identity());
    std::vector<bool> visited(count, false);
    std::queue<int> queue;
    for (int i = 0; i < count; ++i) {
      if (!model.bodies[i].dynamic()) {
        visited[i] = true;
        queue.push(i);
      }
    }
    const auto add = [&](int body, int parent, int edge) {
      JointCoordinates<double> coordinates;
      SE3<double> tree_transform = SE3<double>::Identity();
      if (edge < 0) {
        coordinates = FreeJoint<double>{};
      } else {
        const auto &joint = model.joints[edge];
        const Pose &child_frame =
            joint.a == body ? joint.frame_a : joint.frame_b;
        const Pose &parent_frame =
            joint.a == parent ? joint.frame_a : joint.frame_b;
        frames_[body] = transform(child_frame);
        tree_transform =
            model.bodies[parent].dynamic()
                ? frames_[parent].inverse() * transform(parent_frame)
                : authored_pose(model.initial_state, parent) *
                      transform(parent_frame);
        switch (joint.type) {
        case JointType::Fixed:
          coordinates = FixedJoint<double>{};
          break;
        case JointType::Revolute:
          coordinates = RevoluteJoint<double>{Vec3::UnitZ()};
          break;
        case JointType::Prismatic:
          coordinates = PrismaticJoint<double>{Vec3::UnitZ()};
          break;
        default:
          throw std::invalid_argument("Unsupported ABA joint");
        }
      }
      const auto &source = model.bodies[body];
      const auto &frame = frames_[body];
      const Vec3 com = frame.R.transpose() * (source.com - frame.p);
      const Mat3 inertia = frame.R.transpose() * source.inertia * frame.R;
      const int node = builder.add_body(
          SpatialBody<double>::with_inertia(source.name, source.mass, com,
                                            inertia),
          coordinates, parent < 0 ? 0 : body_to_node_[parent], tree_transform);
      body_to_node_[body] = node;
      node_to_body_.push_back(body);
    };
    const auto drain = [&]() {
      while (!queue.empty()) {
        const int parent = queue.front();
        queue.pop();
        for (int edge : adjacent[parent]) {
          const auto &joint = model.joints[edge];
          const int child = joint.a == parent ? joint.b : joint.a;
          if (visited[child]) {
            continue;
          }
          visited[child] = true;
          add(child, parent, edge);
          queue.push(child);
        }
      }
    };
    drain();
    for (int i = 0; i < count; ++i) {
      if (!visited[i]) {
        visited[i] = true;
        add(i, -1, -1);
        queue.push(i);
        drain();
      }
    }
    tree_ = builder.build();
    set_state(model.initial_state);
  }

  const rapide::State &state() const override { return state_; }
  const StepReport &report() const override { return report_; }

  /** Extract joint configuration and velocity from relative body transforms
   *  and twists, then reconstruct body states with forward kinematics to check
   *  compatibility. No optimization or constraint projection is performed:
   *  inconsistent poses or velocities throw without replacing the old state. */
  void set_state(rapide::State state) override {
    model_.validate_state(state);
    auto coordinates = tree_.zero_state();
    std::vector<SE3<double>> poses(tree_.num_bodies(), SE3<double>::Identity());
    std::vector<Twist<double>> velocities(tree_.num_bodies(),
                                          Twist<double>::Zero());
    for (int node = 1; node < tree_.num_bodies(); ++node) {
      const int body = node_to_body_[node], joint = node - 1;
      const auto &s = state.bodies[body];
      poses[node] = authored_pose(state, body) * frames_[body];
      const Vec3 omega = s.pose.orientation * s.omega;
      const Vec3 arm = s.pose.position - poses[node].p;
      velocities[node] = {poses[node].R.transpose() * omega,
                          poses[node].R.transpose() *
                              (s.velocity - omega.cross(arm))};
      const int parent = tree_.parent[joint];
      const auto relative =
          (poses[parent] * tree_.X_tree[joint]).inverse() * poses[node];
      const auto relative_velocity =
          velocities[node] - (poses[parent].inverse() * poses[node])
                                 .inv_act_on_motion(velocities[parent]);
      const int qi = tree_.q_idx[joint], vi = tree_.v_idx[joint];
      const auto type = joint_type(tree_.joints[joint]);
      if (type == JointType::Free) {
        const Eigen::Quaterniond orientation(relative.R);
        coordinates.q.segment<4>(qi) = orientation.coeffs();
        coordinates.q.segment<3>(qi + 4) = relative.p;
        coordinates.v.segment<3>(vi) = relative_velocity.angular;
        coordinates.v.segment<3>(vi + 3) = relative_velocity.linear;
      } else if (type == JointType::Revolute) {
        coordinates.q(qi) = std::atan2(relative.R(1, 0), relative.R(0, 0));
        coordinates.v(vi) = relative_velocity.angular.z();
      } else if (type == JointType::Prismatic) {
        coordinates.q(qi) = relative.p.z();
        coordinates.v(vi) = relative_velocity.linear.z();
      }
    }
    const auto reconstructed = body_state(coordinates);
    for (int body : node_to_body_) {
      const auto &expected = state.bodies[body];
      const auto &actual = reconstructed.bodies[body];
      if ((expected.pose.position - actual.pose.position).norm() > tolerance_ ||
          expected.pose.orientation.angularDistance(actual.pose.orientation) >
              tolerance_ ||
          (expected.velocity - actual.velocity).norm() > tolerance_ ||
          (expected.omega - actual.omega).norm() > tolerance_) {
        throw std::invalid_argument(
            "ABA initial/body state violates joint constraints");
      }
    }
    coordinates_ = std::move(coordinates);
    state_ = reconstructed;
    report_ = {};
  }

  /** Shift world-frame COM torques to world-origin spatial wrenches for ABA.
   *  Compute generalized accelerations, advance joint velocity/configuration
   *  with shared semi-implicit Euler, and reconstruct public body states.
   *  Tree joints hold by construction, so no constraint projection is needed.
   */
  StepReport advance(double dt, const Forces &forces,
                     const Vec3 &gravity) override {
    std::vector<Wrench<double>> external(tree_.num_bodies(),
                                         Wrench<double>::Zero());
    for (int node = 1; node < tree_.num_bodies(); ++node) {
      const int body = node_to_body_[node];
      const auto &force = forces[body];
      external[node] = {
          force.torque + state_.bodies[body].pose.position.cross(force.force),
          force.force};
    }
    const Eigen::VectorXd effort = Eigen::VectorXd::Zero(tree_.nv);
    const auto acceleration =
        forward_dynamics(tree_, coordinates_, effort, external, gravity);
    auto next = integration::semi_implicit_euler(
        coordinates_.q, coordinates_.v, acceleration, dt,
        [&](const auto &q, const auto &v, double step) {
          return integrate_joint_configuration(tree_, q, v, step);
        });
    JointState<double> coordinates{std::move(next.configuration),
                                   std::move(next.velocity)};
    auto state = body_state(coordinates);
    model_.validate_state(state);
    coordinates_ = std::move(coordinates);
    state_ = std::move(state);
    return report_;
  }

private:
  static SE3<double> transform(const Pose &pose) {
    return {pose.orientation.toRotationMatrix(), pose.position};
  }
  SE3<double> authored_pose(const rapide::State &state, int body) const {
    const auto &s = state.bodies[body];
    const Mat3 rotation = s.pose.orientation.toRotationMatrix();
    return {rotation, s.pose.position - rotation * model_.bodies[body].com};
  }
  /** Map forward kinematics in internal attachment frames back to COM poses,
   *  body-frame angular velocities, and world-frame COM linear velocities.
   *  Static bodies retain their authored state. */
  rapide::State body_state(const JointState<double> &coordinates) const {
    auto result = model_.initial_state;
    const auto kinematics = forward_kinematics(tree_, coordinates);
    for (int node = 1; node < tree_.num_bodies(); ++node) {
      const int body = node_to_body_[node];
      const auto &pose = kinematics.X_world[node];
      const auto &velocity = kinematics.v_body[node];
      const auto &com = tree_.bodies[node].inertia.com;
      auto &s = result.bodies[body];
      s.pose = {pose.p + pose.R * com,
                Eigen::Quaterniond(pose.R * frames_[body].R.transpose())
                    .normalized()};
      s.omega = frames_[body].R * velocity.angular;
      s.velocity = pose.R * (velocity.linear + velocity.angular.cross(com));
    }
    return result;
  }

  rapide::Model model_;
  TreeData<double> tree_;
  JointState<double> coordinates_;
  rapide::State state_;
  StepReport report_;
  std::vector<int> body_to_node_, node_to_body_;
  std::vector<SE3<double>> frames_;
  double tolerance_ = 1e-8;
};

} // namespace rapide::aba
