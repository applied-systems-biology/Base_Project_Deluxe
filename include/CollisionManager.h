//  Copyright by Simone Piccioni
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef COLLISIONMANAGER_H
#define COLLISIONMANAGER_H

#include <random>
#include <geometrycentral/surface/vector_heat_method.h>

#include "../include/Space.h"
#include "../include/Agent.h"


class CollisionManager {
public:
    CollisionManager(Space* s);
    // methods:
    void checkCollisions(std::vector<Agent> agents);
    bool fixCollisions(std::vector<Agent>& agents);


private:
    Space *space;
    std::vector<int> occupation_matrix;

    std::set<std::pair<int, int>> detected_collisions;
    std::map<std::pair<int, int>, std::vector<geometrycentral::Vector2>> resolution_velocities;
    std::map<std::pair<int, int>, double> resolution_lengths;

    std::uniform_real_distribution<double> dist_01;
};



#endif //COLLISIONMANAGER_H
