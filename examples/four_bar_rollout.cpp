#include "examples/models/four_bar.hpp"
#include "rapide/simulator.hpp"
#include "rapide/visualization/recording.hpp"
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>

int main(int argc, char **argv) try {
  using namespace rapide;
  const Pose mount{{0, 0, .8},
                   rotation_increment(Vec3{std::numbers::pi / 2, 0, 0})};
  const auto mechanism = rapide::examples::build_four_bar(.9, mount);
  const auto spring = rapide::examples::four_bar_spring(mechanism);
  Simulator sim(mechanism.model);
  std::unique_ptr<rapide::visualization::Recording> recording;
  if (argc != 1) {
    if (argc != 3 || std::string(argv[1]) != "--record") {
      std::cerr << "Usage: four_bar_rollout [--record recording.json]\n";
      return 1;
    }
    using rapide::visualization::BoxVisual;
    const std::vector<BoxVisual> visuals{
        {mechanism.chassis, {2.2, .14, .12}, {{1, 0, 0}}, 0x566475},
        {mechanism.crank, {1, .075, .075}, {}, 0x3388cc},
        {mechanism.coupler, {2, .075, .075}, {}, 0xe6a23c},
        {mechanism.follower, {1, .075, .075}, {}, 0x36a889}};
    recording = std::make_unique<rapide::visualization::Recording>(
        argv[2], sim.model(), visuals, 50, std::vector<LinearSpring>{spring});
    recording->frame(0, sim.state());
  }
  sim.set_force_callback([spring](const Model &model, const State &state) {
    auto forces = zero_forces(model);
    add_spring(model, state, forces, spring);
    return forces;
  });
  std::cout << "time,angle,angular_velocity,coupler_x,coupler_y,coupler_z,"
               "chassis_z,spring_length,spring_compression,position_error,"
               "orientation_error,velocity_error,rank,"
               "rows\n";
  for (int i = 0; i < 10000; ++i) {
    const auto report = sim.step(.001);
    if (recording && (i + 1) % 20 == 0) {
      recording->frame(sim.time(), sim.state());
    }
    if (i % 10 == 0 || i == 9999) {
      const auto &coupler = sim.state().bodies[mechanism.coupler].pose.position;
      const double length =
          (world_point(sim.model(), sim.state(), spring.b, spring.point_b) -
           world_point(sim.model(), sim.state(), spring.a, spring.point_a))
              .norm();
      std::cout << std::setprecision(10) << sim.time() << ','
                << joint_position(sim.model(), sim.state(), mechanism.drive)
                << ','
                << joint_velocity(sim.model(), sim.state(), mechanism.drive)
                << ',' << coupler.x() << ',' << coupler.y() << ','
                << coupler.z() << ','
                << sim.state().bodies[mechanism.chassis].pose.position.z()
                << ',' << length << ',' << spring.rest_length - length << ','
                << report.position_error << ',' << report.orientation_error
                << ',' << report.velocity_error << ',' << report.rank << ','
                << report.constraint_rows << '\n';
    }
  }
  if (recording) {
    recording->finish();
  }
} catch (const std::exception &error) {
  std::cerr << error.what() << '\n';
  return 1;
}
