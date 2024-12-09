# Building the executable:
-------------------------------------------------------------------------------

CMAKE may be used to build the executable, below is a short shell script that builds a fresh version of the executable.

    #!/bin/bash

    rm -r build
    mkdir build; cd build 
    cmake ../cmake
    PKGS="-D PKG_GRANULAR=on -D PKG_VTK=on -D PKG_SRD=on -D PKG_MPI=on"
    cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake $PKGS ../cmake
    cmake --build . -- -j 10

The GRANULAR package allows DEM simulation in LAMMPS and the VTK package allows creation of VTK files for viewing in paraview. The SRD package is for fluid coupling and the MPI is for parallelization. There are various other ways to build the lammps executable if you are interested: https://docs.lammps.org/Build.html.


# Running the code:
-------------------------------------------------------------------------------

In the sims directory of the code base there is a folder called avicelTableting. Inside there is an input file called in.avicelTableting200. This input file defines a tableting simulation that is composed of die filling, compaction, release, and ejection. The geometry of the simulation is modeled after a real compaction simulator die. The material properties that are set roughly represent Avicel PH102. To run this simulation in serial, navigate inside the avicelTableting folder and run the following command:

    (path to lmp executable) < in.avicelTableting200

The simulation consists of 200 particles and takes ~10 minutes to run on a 2021 M1 Max. To generate a new packing with particles of different radii or distribution you can use the matlab script within generatePackings called generateCylinderPacking.m. It allows you to define the cylinder size, number of desired particles, minimum radius, and maximum radius. It will then generate a packing called spheres.data that is saved in the avicelTableting folder. If you change the distribution of the particle packing make sure to update the following items in the input script:

    neighbor - should be ~1.5x greater than the max radius Rmax.
    timestep - dt = 0.35*sqrt{m/k} = 0.35*sqrt{(rho*4/3*pi*Rmin^3)/(kappa*Rmin)}.
    variable atomRadius - specified value in front of prefactor should match Rmin from generateCylinderPacking.m.

Once running the program should output three files a dump file and two csv files. The dump file can be used in Ovito for visualization. If you desire vtk files for visualization in Paraview create a folder called post and uncomment the dump command in in.avicelTableting. The generated csv files contain the upper punch displacement force relation and particle stress information.

Plots of the axial stress and radial stress versus displacement can be generated using the matlab script avicelTabletingPlotStresses.m. The script will work even while the simulation is in the middle of running and is a good way to check on the status of the simulation.

The other input script in the folder is in.avicelTableting20000 this 20,000 particle simulation is the exact one presented in Zunker et al., 2024 and is most efficiently run in parallel with the following command:

    mpirun --np 4 (path to lmp executable) -in in.avicelTableting20000 

Anticipate this simulation to take longer to run, on a 2021 M1 Max running in parallel on 4 processors the simulation took 28 hours to complete.

Also within the sims directory is the MPFEM folder, which contains two small simulations involving the compaction of 14 monodisperse particles and 12 tridisperse particles as described in Zunker et al., 2024. The input files for these simulations can then be run using the following commands:

    (path to lmp executable) < in.triaxial14particles
    (path to lmp executable) < in.triaxial12particles 

The outputs of each of these simulations are a dump file and a .csv containing the measured forces in each principal direction. To visualize the measured stresses you can run the Matlab script macro_stresses.m. This will plot stresses for both the LAMMPS results and data from corresponding FEM simulations. The Abaqus 2022 input files for the monodisperse and tridisperse FEM are included in the folder as MPFEM.inp and MPFEM_diff_radii.inp, respectively. For convenience, the principal stress output data for these FEM simulations is already included in the folder as fem.mat and fem_diff_radii.mat.   