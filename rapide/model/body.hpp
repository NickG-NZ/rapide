#pragma once
#include <string>
#include "rapide/math/spatial.hpp"

namespace rapide {

template<typename S>
struct Body {
    std::string      name;
    SpatialInertia<S> inertia;

    // Convenience: create a point-mass body.
    static Body point_mass(std::string name, S mass, Eigen::Matrix<S,3,1> com_in_body) {
        return {
            std::move(name),
            {mass, com_in_body, Eigen::Matrix<S,3,3>::Zero()}
        };
    }

    // Create body with full inertia specification.
    static Body with_inertia(std::string name, S mass,
                              Eigen::Matrix<S,3,1> com_in_body,
                              Eigen::Matrix<S,3,3> I_about_com) {
        return {std::move(name), {mass, com_in_body, I_about_com}};
    }
};

}  // namespace rapide
