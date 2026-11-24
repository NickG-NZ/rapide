#include "examples/models/motorcycle.hpp"
#include <iomanip>
#include <iostream>

int main() {
  using namespace rapide;
  using namespace rapide::examples;
  MotorcycleParams<double> parameters;
  const auto bike = build_motorcycle_suspended(parameters);
  Simulator simulator(bike.model, Backend::ABA);
  attach_tire_forces(simulator, bike, parameters);
  attach_suspension(simulator, bike, parameters);
  simulator.set_joint_effort(bike.rear_wheel_joint, 100);
  std::cout << "time,chassis_height,rear_fork,front_fork,forward_velocity\n";
  for (int step = 1; step <= 20000; ++step) {
    simulator.step(.0001);
    if (step % 1000 != 0) {
      continue;
    }
    const auto &chassis = simulator.state().bodies[1];
    std::cout << std::setprecision(8) << simulator.time() << ','
              << chassis.pose.position.z() << ','
              << joint_position(simulator.model(), simulator.state(),
                                bike.rear_fork_joint)
              << ','
              << joint_position(simulator.model(), simulator.state(),
                                bike.front_fork_joint)
              << ',' << chassis.velocity.x() << '\n';
  }
}
