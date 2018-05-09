Grupo de Modelacion Matematica y Computacional
Instituto de Goefisica
Universidad Nacional Autonoma de Mexico
Mexico D.F.

Autores:
Ismael Herrera Revilla
Robert Yates Smit
Antonio Carrillo Ledesma

Codigo liberado bajo la licencia GPL ver. 2.0


Si se quiere usar make
  $ ./Link
  $ make deps
  $ make
correr con make de ejemplo por omsion  
  $ make run
correr usando MPICH 1.x
  $ mpirun.mpich -np 4 DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
correr usando MPICH 2.x
  $ mpirun.mpich2 -np 4 DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
correr usando MPI
  $ mpirun -np 4 DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10
Limpiar archivos de compilación
  $ make clean
  $ ./UnLink



Si no se quiere usar make, para compilar los fuentes respectivos
  $ ./Link
compilar usando MPICH 1.x
  $ mpiCC.mpich -O2 *.cpp -o DVS-DDM -lm
compilar usando MPICH 2.x
  $ mpic++.mpich2 -O2 *.cpp -o DVS-DDM -lm
compilar usando MPI
  $ mpic++ -O2 *.cpp -o DVS-DDM -lm

ejecutar en Linux:
  $ lamboot -v
correr usando MPICH 1.x
  $ mpirun.mpich -np 4 DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
correr usando MPICH 2.x
  $ mpirun.mpich2 -np 4 DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
correr usando MPI
  $ mpirun -np 4 DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10
  $ lamhalt -v
  $ rm DVS-DDM
  $ ./Unlink

