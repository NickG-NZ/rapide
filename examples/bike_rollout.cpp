#include "examples/models/bike.hpp"
#include "rapide/simulator.hpp"
#include "rapide/visualization/recording.hpp"
#include <iomanip>
#include <iostream>
#include <memory>

int main(int argc, char **argv) try {
  using namespace rapide;
  using namespace rapide::visualization;
  const auto bike = rapide::examples::build_bike(2, 0, .15);
  Simulator sim(bike.model);
  std::unique_ptr<Recording> recording;
  if (argc != 1) {
    if (argc != 3 || std::string(argv[1]) != "--record") {
      std::cerr << "Usage: bike_rollout [--record recording.json]\n";
      return 1;
    }
    std::vector<BoxVisual> visuals{{0, {14, 6, .04}, {{5, 0, -.02}}, 0x718096}};
    const auto bar = [&](BodyId body, const Vec3 &a, const Vec3 &b,
                         unsigned color) {
      visuals.push_back({body,
                         {(b - a).norm(), .045, .045},
                         {(a + b) / 2, Eigen::Quaterniond::FromTwoVectors(
                                           Vec3::UnitX(), b - a)},
                         color});
    };
    const Vec3 pivot{-.2, 0, -.17}, seat{-.15, 0, .12}, head{.4, 0, .12},
        bottom{-.1, 0, -.28};
    const Vec3 front = bike.model.joints.back().frame_a.position;
    for (const auto &edge :
         std::vector<std::pair<Vec3, Vec3>>{{pivot, seat},
                                            {pivot, bottom},
                                            {seat, bottom},
                                            {seat, head},
                                            {bottom, head},
                                            {head, front},
                                            {bottom, bike.spring.point_a}}) {
      bar(bike.chassis, edge.first, edge.second, 0xE68A32);
    }
    bar(bike.swingarm, {-.45, 0, -.15}, {0, 0, 0}, 0x3388CC);
    bar(bike.swingarm, {0, 0, 0}, bike.spring.point_b, 0x3388CC);
    visuals.push_back(
        {bike.chassis, {.28, .18, .05}, {{-.15, 0, .16}}, 0x23364A});
    bar(bike.chassis, head, {.4, 0, .25}, 0xAAB7C4);
    bar(bike.chassis, {.4, -.25, .25}, {.4, .25, .25}, 0x23364A);
    std::vector<TireVisual> tires;
    for (auto wheel : {bike.rear, bike.front}) {
      tires.push_back({wheel, bike.tire.radius, bike.tire.rim_radius,
                       bike.tire.width, bike.tire.ground_height,
                       bike.tire.ground_normal});
    }
    recording = std::make_unique<Recording>(
        argv[2], sim.model(), visuals, 50,
        std::vector<LinearSpring>{bike.spring}, tires);
    recording->frame(0, sim.state());
  }
  sim.set_force_callback([&](const Model &model, const State &state) {
    return rapide::examples::bike_forces(bike, model, state);
  });
  std::cout << "time,x,y,z,vx,vy,lean,rear_lateral_slip,front_lateral_slip,"
               "rear_lateral_force,front_lateral_force,position_error,"
               "suspension_angle,spring_length,spring_compression\n";
  for (int step = 1; step <= 3000; ++step) {
    const auto report = sim.step(.001);
    if (step % 20 != 0) {
      continue;
    }
    if (recording) {
      recording->frame(sim.time(), sim.state());
    }
    const auto &frame = sim.state().bodies[bike.chassis];
    const auto rear = bike.tire.compute(sim.state().bodies[bike.rear]);
    const auto front = bike.tire.compute(sim.state().bodies[bike.front]);
    const Vec3 axle = frame.pose.orientation * Vec3::UnitY();
    const double lean = std::asin(std::clamp(axle.z(), -1., 1.));
    const double spring_length =
        (world_point(sim.model(), sim.state(), bike.spring.b,
                     bike.spring.point_b) -
         world_point(sim.model(), sim.state(), bike.spring.a,
                     bike.spring.point_a))
            .norm();
    std::cout << std::setprecision(10) << sim.time() << ','
              << frame.pose.position.x() << ',' << frame.pose.position.y()
              << ',' << frame.pose.position.z() << ',' << frame.velocity.x()
              << ',' << frame.velocity.y() << ',' << lean << ','
              << rear.lateral_slip << ',' << front.lateral_slip << ','
              << rear.lateral_force << ',' << front.lateral_force << ','
              << report.position_error << ','
              << joint_position(sim.model(), sim.state(), bike.suspension)
              << ',' << spring_length << ','
              << bike.spring.rest_length - spring_length << '\n';
    if (std::abs(lean) > .7) {
      std::cerr << "Stopped at " << sim.time()
                << " s: lean exceeded 0.7 rad; frame-ground contact is not "
                   "modeled.\n";
      break;
    }
  }
  if (recording) {
    recording->finish();
  }
} catch (const std::exception &error) {
  std::cerr << error.what() << '\n';
  return 1;
}
