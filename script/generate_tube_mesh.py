#  Copyright by Simone Piccioni
#  Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
#  https://www.leibniz-hki.de/en/applied-systems-biology.html
#  HKI-Center for Systems Biology of Infection
#  Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
#  Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
#  This code is licensed under BSD 2-Clause
#  See the LICENSE file provided with this code for the full license.

import gmsh

tube_length = 10
tube_radius = 2

gmsh.initialize()
gmsh.model.add("tube")

tube = gmsh.model.occ.addCylinder(0,0,0,
                                  0,0,tube_length,
                                  tube_radius)
gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)

gmsh.model.mesh.generate(2)

gmsh.write("../input/mesh/tube_mesh.stl")

gmsh.fltk.run()
gmsh.finalize()


