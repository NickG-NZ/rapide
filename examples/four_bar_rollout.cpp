#include "examples/models/four_bar.hpp"
#include "rapide/maximal/simulator.hpp"
#include <iomanip>
#include <iostream>

int main() {
  using namespace rapide::maximal;
  const auto mechanism = rapide::examples::build_four_bar();
  Simulator sim(mechanism.model);
  sim.gravity.setZero();
  sim.set_force_callback(
      [drive = mechanism.drive](const Model &model, const State &state) {
        auto forces = zero_forces(model);
        add_joint_effort(model, state, forces, drive,
                         .5 - .2 * joint_velocity(model, state, drive));
        return forces;
      });
  std::cout << "time,angle,position_error,orientation_error,velocity_error,"
               "rank,rows\n";
  for (int i = 0; i < 5000; ++i) {
    const auto report = sim.step(.001);
    if (i % 100 == 0)
      std::cout << std::setprecision(10) << sim.time() << ','
                << joint_position(sim.model(), sim.state(), mechanism.drive)
                << ',' << report.position_error << ','
                << report.orientation_error << ',' << report.velocity_error
                << ',' << report.rank << ',' << report.constraint_rows << '\n';
  }
}
