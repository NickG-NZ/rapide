#pragma once
#include "rapide/math/pose.hpp"
#include "rapide/model/forces.hpp"
#include <fstream>
#include <iomanip>

namespace rapide::visualization {

struct BoxVisual {
  rapide::BodyId body;
  rapide::Vec3 size;
  rapide::Pose local;
  unsigned color = 0x4a90d9;
};

/** Elliptic torus tire centered at the body's COM, with its axle along local Y.
 */
struct TireVisual {
  rapide::BodyId body;
  double radius, rim_radius, width, ground_height = 0;
  rapide::Vec3 ground_normal = rapide::Vec3::UnitZ();
};

/** Stream a versioned scene and COM poses for viewer-independent playback.
 *  Call finish() after the final sample; quaternions are serialized as xyzw. */
class Recording {
public:
  Recording(const std::string &path, const rapide::Model &model,
            const std::vector<BoxVisual> &visuals, double fps,
            const std::vector<rapide::LinearSpring> &springs = {},
            const std::vector<TireVisual> &tires = {})
      : output_(path), body_count_(model.bodies.size()) {
    if (!output_) {
      throw std::runtime_error("Cannot open recording: " + path);
    }
    if (!std::isfinite(fps) || fps <= 0) {
      throw std::invalid_argument("Invalid recording rate");
    }
    output_.exceptions(std::ios::badbit | std::ios::failbit);
    output_ << std::setprecision(17) << "{\"version\":1,\"fps\":" << fps
            << ",\"bodies\":[";
    for (size_t i = 0; i < model.bodies.size(); ++i) {
      if (i) {
        output_ << ',';
      }
      output_ << "{\"name\":";
      string(model.bodies[i].name);
      output_ << ",\"com\":";
      vector(model.bodies[i].com);
      output_ << '}';
    }
    output_ << "],\"visuals\":[";
    for (size_t i = 0; i < visuals.size(); ++i) {
      const auto &v = visuals[i];
      if (v.body < 0 || static_cast<size_t>(v.body) >= body_count_ ||
          !v.size.allFinite() || (v.size.array() <= 0).any()) {
        throw std::invalid_argument("Invalid visual geometry");
      }
      rapide::validate_pose(v.local);
      if (i) {
        output_ << ',';
      }
      output_ << "{\"body\":" << v.body << ",\"type\":\"box\",\"size\":";
      vector(v.size);
      output_ << ",\"pose\":";
      pose(v.local);
      output_ << ",\"color\":" << v.color << '}';
    }
    output_ << "],\"joints\":[";
    for (size_t i = 0; i < model.joints.size(); ++i) {
      const auto &j = model.joints[i];
      if (i) {
        output_ << ',';
      }
      output_ << "{\"a\":" << j.a << ",\"b\":" << j.b << ",\"point_a\":";
      vector(j.frame_a.position);
      output_ << ",\"point_b\":";
      vector(j.frame_b.position);
      output_ << '}';
    }
    output_ << "],\"springs\":[";
    for (size_t i = 0; i < springs.size(); ++i) {
      const auto &s = springs[i];
      if (s.a < 0 || s.b < 0 || static_cast<size_t>(s.a) >= body_count_ ||
          static_cast<size_t>(s.b) >= body_count_ || !s.point_a.allFinite() ||
          !s.point_b.allFinite() || !std::isfinite(s.stiffness) ||
          s.stiffness < 0 || !std::isfinite(s.damping) || s.damping < 0 ||
          !std::isfinite(s.rest_length) || s.rest_length <= 0) {
        throw std::invalid_argument("Invalid spring visual");
      }
      if (i) {
        output_ << ',';
      }
      output_ << "{\"a\":" << s.a << ",\"b\":" << s.b << ",\"point_a\":";
      vector(s.point_a);
      output_ << ",\"point_b\":";
      vector(s.point_b);
      output_ << ",\"stiffness\":" << s.stiffness
              << ",\"damping\":" << s.damping
              << ",\"rest_length\":" << s.rest_length << '}';
    }
    output_ << "],\"tires\":[";
    for (size_t i = 0; i < tires.size(); ++i) {
      const auto &t = tires[i];
      if (t.body < 0 || static_cast<size_t>(t.body) >= body_count_ ||
          !std::isfinite(t.radius) || !std::isfinite(t.rim_radius) ||
          !std::isfinite(t.width) || !std::isfinite(t.ground_height) ||
          !t.ground_normal.allFinite() ||
          std::abs(t.ground_normal.norm() - 1) > 1e-8 || t.rim_radius <= 0 ||
          t.radius <= t.rim_radius || t.width <= 0) {
        throw std::invalid_argument("Invalid tire visual");
      }
      if (i) {
        output_ << ',';
      }
      output_ << "{\"body\":" << t.body << ",\"radius\":" << t.radius
              << ",\"rim_radius\":" << t.rim_radius << ",\"width\":" << t.width
              << ",\"ground_height\":" << t.ground_height
              << ",\"ground_normal\":";
      vector(t.ground_normal);
      output_ << '}';
    }
    output_ << "],\"frames\":[";
  }

  void frame(double time, const rapide::State &state) {
    if (finished_ || !std::isfinite(time) || time < 0 ||
        (!first_ && time <= last_time_) || state.bodies.size() != body_count_) {
      throw std::invalid_argument("Invalid recording frame");
    }
    for (const auto &body : state.bodies) {
      rapide::validate_pose(body.pose);
    }
    if (!first_) {
      output_ << ',';
    }
    first_ = false;
    last_time_ = time;
    output_ << "{\"time\":" << time << ",\"poses\":[";
    for (size_t i = 0; i < state.bodies.size(); ++i) {
      if (i) {
        output_ << ',';
      }
      pose(state.bodies[i].pose);
    }
    output_ << "]}";
  }

  void finish() {
    if (!finished_) {
      output_ << "]}\n";
      output_.close();
      finished_ = true;
    }
  }

private:
  void string(const std::string &value) {
    output_ << '"';
    for (unsigned char c : value) {
      if (c == '"' || c == '\\') {
        output_ << '\\' << c;
      } else if (c < 32) {
        output_ << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                << unsigned(c) << std::dec;
      } else {
        output_ << c;
      }
    }
    output_ << '"';
  }
  void vector(const rapide::Vec3 &v) {
    output_ << '[' << v.x() << ',' << v.y() << ',' << v.z() << ']';
  }
  void pose(const rapide::Pose &p) {
    output_ << '[' << p.position.x() << ',' << p.position.y() << ','
            << p.position.z() << ',' << p.orientation.x() << ','
            << p.orientation.y() << ',' << p.orientation.z() << ','
            << p.orientation.w() << ']';
  }
  std::ofstream output_;
  size_t body_count_;
  bool first_ = true, finished_ = false;
  double last_time_ = 0;
};

} // namespace rapide::visualization
