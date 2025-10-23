# BASE PROJECT DELUXE
Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
https://www.leibniz-hki.de/en/applied-systems-biology.html
HKI-Center for Systems Biology of Infection
Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
Adolf-Reichwein-Straße 23, 07745 Jena, Germany

This code is licensed under BSD 2-Clause
See the LICENSE file provided with this code for the full license.

---

This project is a basic example of an hABM (agents + PDEs) simulation, where "good" cells scan a 2D-surface (defined by a planar mesh) using BPRW (see [link](https://varennes.github.io/bprw/)), in search of "evil" cells. The good cells are guided by chemical signals (chemokines) released by evil cells.

Example of a simulation on a planar surface:

![alooo](imgs/simulation_plane.gif)

## LIBRARIES & DEPENDENCIES
This project uses the following libraries:

- [geometry-central](https://geometry-central.net/) for intrinsic geodesic agent motion and surface metrics
- [gmsh](https://gmsh.info/) to create meshes (using the Python wrapper)
- [MFEM](https://mfem.org/) for solving PDEs
- (optional) [ParaView](https://www.paraview.org/) for visualizing simulation outputs

## INSTALLATION & RUNNING
(NOTE: The following steps are meant for Ubuntu 24, with an other OS you might have to build the libraries from the source)

Basic requirements:
```bash
sudo apt install git g++ cmake
```

Install gmsh:
```bash
sudo apt install gmsh libgmsh-dev
```

Install MFEM:
```bash
sudo apt install libmfem-dev
```

Install geometry-central dependencies:
```bash
sudo apt install libeigen3-dev libsuitesparse-dev
```

Install geometry-central from source:
```bash
git clone --recurse-submodules https://github.com/nmwsharp/geometry-central.git && cd geometry-central
mkdir build && cd build
cmake ..
make -j4
sudo make install #to install the library in /usr/local/include
```

Compile Base_Project_Deluxe:
```bash
git clone https://github.com/PiSimo/Base_Project_Deluxe.git && cd Base_Project_Deluxe
mkdir build && cd build
cmake ..
make -j4
```

Running:
```bash
./Base_Project_Deluxe
```
the output files will be saved in <code>output/</code>.
## OVERVIEW
![overview](imgs/project_overview.jpg)

more details [here](project_overview.ipynb).
