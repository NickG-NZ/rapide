#pragma once
#include "rapide/maximal/state.hpp"
#include "rapide/model/model.hpp"
#include <algorithm>

namespace rapide::maximal {

/** Solver row layout for a physical joint; independent of model topology. */
struct JointConstraints {
  int linear_rows;
  int angular_rows;

  explicit JointConstraints(JointType type) {
    switch (type) {
    case JointType::Ball:
      linear_rows = 3;
      angular_rows = 0;
      break;
    case JointType::Fixed:
      linear_rows = 3;
      angular_rows = 3;
      break;
    case JointType::Revolute:
      linear_rows = 3;
      angular_rows = 2;
      break;
    case JointType::Prismatic:
      linear_rows = 2;
      angular_rows = 3;
      break;
    case JointType::Free:
      throw std::invalid_argument("Maximal free motion needs no joint");
    default:
      throw std::invalid_argument("Unknown joint type");
    }
  }

  int rows() const { return linear_rows + angular_rows; }
};

struct Constraints {
  Eigen::VectorXd residual;
  Eigen::MatrixXd jacobian;
  std::vector<bool> angular;

  double position_error() const { return error(false); }
  double orientation_error() const { return error(true); }

private:
  double error(bool rotation) const {
    double maximum = 0;
    for (int i = 0; i < residual.size(); ++i) {
      if (angular[i] == rotation) {
        maximum = std::max(maximum, std::abs(residual(i)));
      }
    }
    return maximum;
  }
};

/** Assemble phi and J such that phi_dot = J [omega_body; v_com_world].
 *  Orientation residuals use the principal SO(3) logarithm; its chart is
 *  discontinuous at a relative rotation of pi. */
inline Constraints assemble_constraints(const Model &model,
                                        const State &state) {
  int rows = 0;
  for (const auto &joint : model.joints) {
    rows += JointConstraints(joint.type).rows();
  }
  Constraints out{Eigen::VectorXd::Zero(rows),
                  Eigen::MatrixXd::Zero(rows, model.nv),
                  std::vector<bool>(rows, false)};
  int row = 0;
  for (const auto &joint : model.joints) {
    const JointConstraints layout(joint.type);
    const auto &a = state.bodies[joint.a];
    const auto &b = state.bodies[joint.b];
    const Mat3 ra = a.pose.orientation.toRotationMatrix();
    const Mat3 rb = b.pose.orientation.toRotationMatrix();
    const Vec3 arm_a =
        ra * (joint.frame_a.position - model.bodies[joint.a].com);
    const Vec3 arm_b =
        rb * (joint.frame_b.position - model.bodies[joint.b].com);
    const Vec3 d = b.pose.position + arm_b - a.pose.position - arm_a;
    const Eigen::Quaterniond qa =
        a.pose.orientation * joint.frame_a.orientation;
    const Eigen::Quaterniond qb =
        b.pose.orientation * joint.frame_b.orientation;
    Eigen::Matrix<double, 6, 6> ja = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 6, 6> jb = ja;
    ja.block<3, 3>(0, 0) = skew(arm_a) * ra;
    ja.block<3, 3>(0, 3) = -Mat3::Identity();
    jb.block<3, 3>(0, 0) = -skew(arm_b) * rb;
    jb.block<3, 3>(0, 3) = Mat3::Identity();
    const int linear_rows = layout.linear_rows;
    if (joint.type == JointType::Prismatic) {
      const Mat3 t = qa.toRotationMatrix().transpose();
      ja.block<3, 3>(0, 0) = (t * (ja.block<3, 3>(0, 0) + skew(d) * ra)).eval();
      ja.block<3, 3>(0, 3) = -t;
      jb.topRows<3>() = (t * jb.topRows<3>()).eval();
      out.residual.segment<2>(row) = (t * d).head<2>();
    } else {
      out.residual.segment<3>(row) = d;
    }
    if (layout.angular_rows != 0) {
      const Vec3 e = rotation_log(qa.conjugate() * qb);
      const Mat3 t =
          inverse_left_jacobian(e) * qa.toRotationMatrix().transpose();
      const int nr = layout.angular_rows;
      ja.block(linear_rows, 3, nr, 3).setZero();
      jb.block(linear_rows, 3, nr, 3).setZero();
      ja.block(linear_rows, 0, nr, 3) = (-t * ra).topRows(nr);
      jb.block(linear_rows, 0, nr, 3) = (t * rb).topRows(nr);
      out.residual.segment(row + linear_rows, nr) = e.head(nr);
      for (int k = 0; k < nr; ++k) {
        out.angular[row + linear_rows + k] = true;
      }
    }
    const int oa = model.bodies[joint.a].velocity_offset;
    const int ob = model.bodies[joint.b].velocity_offset;
    if (oa >= 0) {
      out.jacobian.block(row, oa, layout.rows(), 6) = ja.topRows(layout.rows());
    }
    if (ob >= 0) {
      out.jacobian.block(row, ob, layout.rows(), 6) = jb.topRows(layout.rows());
    }
    row += layout.rows();
  }
  return out;
}

} // namespace rapide::maximal
