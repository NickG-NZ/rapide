#include "examples/models/planar_wheel.hpp"
#include "rapide/simulator.hpp"
#include "rapide/visualization/recording.hpp"
#include <iomanip>
#include <iostream>
#include <memory>

int main(int argc, char **argv) try {
  using namespace rapide;
  using namespace rapide::visualization;
  std::string recording_path;
  double camber = 0;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--record" && i + 1 < argc) {
      recording_path = argv[++i];
    } else if (arg == "--camber" && i + 1 < argc) {
      const std::string value = argv[++i];
      size_t consumed = 0;
      camber = std::stod(value, &consumed);
      if (consumed != value.size() || !std::isfinite(camber)) {
        throw std::invalid_argument("Camber must be a finite angle in radians");
      }
    } else {
      std::cerr << "Usage: wheel_rollout [--record recording.json] [--camber "
                   "radians]\n";
      return 1;
    }
  }
  const auto mechanism = rapide::examples::build_planar_wheel(.47, .65, camber);
  Simulator sim(mechanism.model);
  std::unique_ptr<Recording> recording;
  if (!recording_path.empty()) {
    const std::vector<BoxVisual> visuals{
        {0, {6, 1.2, .04}, {{1.5, 0, -.02}}, 0x718096}};
    recording = std::make_unique<Recording>(
        recording_path, sim.model(), visuals, 50, std::vector<LinearSpring>{},
        std::vector<TireVisual>{
            {mechanism.wheel, mechanism.tire.radius, mechanism.tire.rim_radius,
             mechanism.tire.width, mechanism.tire.ground_height,
             mechanism.tire.ground_normal}});
    recording->frame(0, sim.state());
  }
  sim.set_force_callback([&](const Model &model, const State &state) {
    auto forces = zero_forces(model);
    forces[mechanism.wheel] =
        mechanism.tire.compute(state.bodies[mechanism.wheel]).wrench;
    return forces;
  });
  std::cout << "time,x,height,vx,vz,omega,compression,slip,normal_force,"
               "force_x,patch_x,patch_y,patch_z,longitudinal_slip,lateral_slip,"
               "longitudinal_force,lateral_force\n";
  for (int i = 1; i <= 6000; ++i) {
    sim.step(.001);
    if (recording && i % 20 == 0) {
      recording->frame(sim.time(), sim.state());
    }
    if (i % 10 == 0) {
      const auto &s = sim.state().bodies[mechanism.wheel];
      const auto tire = mechanism.tire.compute(s);
      std::cout << std::setprecision(10) << sim.time() << ','
                << s.pose.position.x() << ',' << s.pose.position.z() << ','
                << s.velocity.x() << ',' << s.velocity.z() << ',' << s.omega.y()
                << ',' << tire.compression << ',' << tire.slip << ','
                << tire.wrench.force.z() << ',' << tire.wrench.force.x() << ','
                << tire.contact_point.x() << ',' << tire.contact_point.y()
                << ',' << tire.contact_point.z() << ','
                << tire.longitudinal_slip << ',' << tire.lateral_slip << ','
                << tire.longitudinal_force << ',' << tire.lateral_force << '\n';
    }
  }
  if (recording) {
    recording->finish();
  }
} catch (const std::exception &error) {
  std::cerr << error.what() << '\n';
  return 1;
}
