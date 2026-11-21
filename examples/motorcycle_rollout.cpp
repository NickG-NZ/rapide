// Demonstrate: motorcycle with tire contact + suspension, driven forward.
// Prints chassis height, fork compression, and forward speed each 100ms.
#include <iostream>
#include <iomanip>
#include "rapide/simulator.hpp"
#include "examples/models/motorcycle.hpp"

using namespace rapide;
using namespace rapide::examples;
using S = double;

int main() {
    MotorcycleParams<S> p;
    p.has_suspension = true;

    // ── Build model with suspension ───────────────────────────────────────────
    auto mm    = build_motorcycle_suspended(p);
    auto& model = mm.model;

    std::cout << "Bodies: " << model.num_bodies() << "\n";
    std::cout << "nq=" << model.nq << "  nv=" << model.nv << "\n\n";

    // Joint/body layout:
    //   0: world        —  FreeJoint(0)     → chassis(1)
    //   1: chassis      —  PrismaticJoint(1) → rear_slider(2)
    //   2: rear_slider  —  RevoluteJoint(2)  → rear_wheel(3)
    //   3: rear_wheel   —  PrismaticJoint(3) → front_slider(4)
    //   4: front_slider —  RevoluteJoint(4)  → front_wheel(5)
    //   5: front_wheel
    for (int i = 0; i < model.num_bodies(); ++i)
        std::cout << "  body[" << i << "] " << model.bodies[i].name << "\n";
    std::cout << "\n";

    // ── Build simulator ───────────────────────────────────────────────────────
    Simulator<S> sim(model);  // copy model into sim
    sim.set_state(motorcycle_standing_state(sim.model(), p.wheel_radius));

    // Tire forces
    attach_tire_forces(sim, mm, p);

    // Suspension springs
    attach_suspension(sim, mm, p);

    // Drive torque on rear wheel (joint index 2 → rear_wheel revolute)
    // v_idx for rear_wheel revolute = v_idx[2]
    const int rear_wheel_joint_idx = 2;  // joint that drives rear wheel
    const S drive_torque = S(100);       // [N·m]

    // ── Header ────────────────────────────────────────────────────────────────
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  t[s]  chassis_z[m]  rear_fork[m]  front_fork[m]  vx[m/s]\n";
    std::cout << "  ----  ------------  ------------  -------------  -------\n";

    const S dt   = S(1e-4);    // 10 kHz — stiff tire spring needs small dt
    const S t_end = S(2.0);
    const int print_every = 1000;  // print at 10Hz

    for (int step = 0; sim.time < t_end; ++step) {
        // Set drive torque
        Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(sim.model().nv);
        tau(sim.model().v_idx[rear_wheel_joint_idx]) = drive_torque;
        sim.set_tau(tau);

        sim.step(dt);

        if (step % print_every == 0) {
            const State<S>& st = sim.state();

            // Chassis position: q[4..6] = [qx,qy,qz,qw,px,py,pz] → pz = q(6)
            S chassis_z = st.q(6);

            // Fork compression: q_fork at joint idx 1 (rear) and 3 (front)
            S rear_fork  = st.q(sim.model().q_idx[mm.rear_fork_joint]);
            S front_fork = st.q(sim.model().q_idx[mm.front_fork_joint]);

            // Forward speed: chassis v_x = v(3) (free joint v slice: [wx,wy,wz,vx,vy,vz])
            S vx = st.v(3);  // v_idx[free_joint]=0, vx at offset 3

            std::cout << "  " << std::setw(6) << sim.time
                      << "  " << std::setw(12) << chassis_z
                      << "  " << std::setw(12) << rear_fork
                      << "  " << std::setw(13) << front_fork
                      << "  " << std::setw(7)  << vx << "\n";
        }
    }

    return 0;
}
