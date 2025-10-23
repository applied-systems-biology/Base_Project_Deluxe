//  Copyright by Simone Piccioni
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#include "../include/CollisionManager.h"
#include <geometrycentral/surface/exact_geodesics.h>

#define REMOVAL_PROBABILITY .7

using namespace std;

CollisionManager::CollisionManager(Space* s): space(s), occupation_matrix(s->gc_mesh->nFaces(),0), dist_01(0.0, 1.0){}

void CollisionManager::checkCollisions(std::vector<Agent> agents) {
    // re-set occupation matrix:
    std::fill(occupation_matrix.begin(), occupation_matrix.end(), 0);

    // remove passed events:
    detected_collisions.clear();
    resolution_velocities.clear();
    resolution_lengths.clear();
    
    // cycle over the agents:
    for(int n_agent_a=0; n_agent_a <agents.size(); n_agent_a++) {
        unordered_set<geometrycentral::surface::Face> faces = agents[n_agent_a].findFacesWithinRadius();


        // now for each face we check if it's occupied by more than 1 agent:
        for (auto f: faces) {
            int face_index = f.getIndex();
            if (occupation_matrix[face_index] == 0 || agents[n_agent_a].getAgentId() == occupation_matrix[face_index]) {
                occupation_matrix[face_index] = agents[n_agent_a].getAgentId();
            }
            else {
                // find the index of the other agent inside of "agents"
                int n_agent_b = -1;
                for (int r=0;r<agents.size();r++) {
                    if (agents[r].getAgentId() == occupation_matrix[face_index]) {
                        n_agent_b = r;
                        break;
                    }
                }
                // NOTE: it's reversed because n_agent_b should be < n_agent_a
                auto candidate_indices = std::make_pair(n_agent_b, n_agent_a);


                // if we didn't already resolve the agent-agent collision:
                if (!detected_collisions.contains(candidate_indices)) {
                    detected_collisions.insert(candidate_indices);

                    geometrycentral::surface::SurfacePoint start = agents[n_agent_a].gc_position;
                    geometrycentral::surface::SurfacePoint end = agents[n_agent_b].gc_position;

                    double radii_sum = agents[n_agent_a].getAgentRadius() + agents[n_agent_b].getAgentRadius();

                    // local MMP algorithm:
                    geometrycentral::surface::GeodesicAlgorithmExact geodesic_solver(*space->gc_mesh, *space->gc_geometry);
                    geodesic_solver.propagate(start,
                                                radii_sum, {end});

                    double total_length = geodesic_solver.getDistance(end);
                    double length_difference = (agents[n_agent_b].getAgentRadius() + agents[n_agent_a].getAgentRadius()) - total_length;

                    if (length_difference < 0) {
                        // the two agents are not touching!
                        resolution_lengths[candidate_indices] = 0;
                        continue;
                    }
                    resolution_lengths[candidate_indices] = length_difference;

                    // the two agents are touching, compute the geodesic to get the resolving directions:

                    geometrycentral::Vector2 resolving_direction_agent_a = (-1) *geodesic_solver.getLog(end).normalize();
                    geometrycentral::Vector2 resolving_direction_agent_b = geodesic_solver.getDistanceGradient(end).normalize(); // this will already be in agent b's tangent plane

                    resolution_velocities[candidate_indices] = {resolving_direction_agent_b, resolving_direction_agent_a};
                }
            }
        }
    }
}

// is true if some agents have been removed
bool CollisionManager::fixCollisions(std::vector<Agent>& agents) {
    std::unordered_set<int> agent_index_killing_list;

    // resolve collisions:
    for (auto collision : detected_collisions) {
        double corrective_displacement = resolution_lengths[{collision.first, collision.second}];

        auto& agent_1 = agents[collision.first];
        auto& agent_2 = agents[collision.second];

        // valid good/good (or evil/evil) interaction: bouncing off
        if (corrective_displacement > 0.0 && agent_1.getAgentType()*agent_2.getAgentType() > 0) {
            geometrycentral::Vector2 old_dir_1 = agent_1.gc_local_velocity_direction;
            geometrycentral::Vector2 old_dir_2 = agent_2.gc_local_velocity_direction;

            agent_1.gc_local_velocity_direction = resolution_velocities[collision][0].normalize();
            agent_2.gc_local_velocity_direction = resolution_velocities[collision][1].normalize();

            double angle_1 = geometrycentral::angle(old_dir_1, agent_1.gc_local_velocity_direction);
            double angle_2 = geometrycentral::angle(old_dir_2, agent_2.gc_local_velocity_direction);

            // the "move" also updates the local velocity in the new reference frame:
            agent_1.move(corrective_displacement/2,false);
            agent_2.move(corrective_displacement/2, false);



            //if you remove the following lines you have elastic collisions:

            // now we rotate back to the original velocity (this is now in the new ref.)
            agent_1.gc_local_velocity_direction = agent_1.gc_local_velocity_direction.rotate(angle_1);
            agent_2.gc_local_velocity_direction = agent_2.gc_local_velocity_direction.rotate(angle_2);

        }
        // valid good/evil (=evil/good) interaction: possible phagocytosis
        else if (corrective_displacement > 0.0 && agent_1.getAgentType()*agent_2.getAgentType()<0) {
            geometrycentral::Vector2 old_dir_1 = agent_1.gc_local_velocity_direction;
            geometrycentral::Vector2 old_dir_2 = agent_2.gc_local_velocity_direction;

            // make the good cell go even closer:
            //if (agent_1.getAgentType() > 0) agent_1.gc_local_velocity_direction = -resolution_velocities[collision][0].normalize();
            //else agent_2.gc_local_velocity_direction = -resolution_velocities[collision][1].normalize();

            // phagocytose:
            if (float p = dist_01(space->gen); p < REMOVAL_PROBABILITY) {
                if (agent_1.getAgentType() < 0) {
                    agent_index_killing_list.insert(agent_1.getAgentId());
                }else agent_index_killing_list.insert(agent_2.getAgentId());
            }
        }
    }

    if (agent_index_killing_list.size() >0) {
        agents.erase(
            std::remove_if(
                agents.begin(),
                agents.end(),
                [&]( Agent& a) {
                    // condition: should this element be removed?
                    return agent_index_killing_list.contains(a.getAgentId());
                }
            ),
            agents.end()
        );
        return true;
    }
    return false;
}
