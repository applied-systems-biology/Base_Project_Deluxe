//  Copyright by Simone Piccioni
//  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
//  https://www.leibniz-hki.de/en/applied-systems-biology.html
//  HKI-Center for Systems Biology of Infection
//  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
//  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
//  This code is licensed under BSD 2-Clause
//  See the LICENSE file provided with this code for the full license.

#ifndef UTILS_H
#define UTILS_H

#include "Space.h"
#include "Agent.h"
#include "Field.h"




namespace utils {
    void saveAgentsPositionToFile(Space* space,std::vector<Agent>& agents, std::string file_name);
    void saveAgentsFacesToFile(Space* space, std::vector<Agent>& agents, std::string file_name);

    void saveFieldToFile(Space* space, Field& field, std::string file_name, std::string field_name="scalar field");


}
#endif //UTILS_H
