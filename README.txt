Métodos de Descomposición de Dominio en el Espacio De Vectores Derivados y su Implementación Computacional en Paralelo 

Grupo de Modelacion Matematica y Computacional
Instituto de Goefisica
Universidad Nacional Autonoma de Mexico
Mexico D.F.

Autores:
 Ismael Herrera Revilla
 Robert Yates Smit
 Antonio Carrillo Ledesma
 Alberto Rosas Medina


Codigo liberado bajo la licencia GPL ver. 2.0



Fuentes del codigo del metodo de DVS-DDM secuencial como paralelo, el cual se separo en una estructura
de directorios para un manejo mas facil, la estructura queda de la siguiente manera:


DVS                        Implementacion de DVS-DDM
Methods                    Implementacion de los distintos metodos de DVS
Library                    Biblioteca de funciones
Problems                   Implementacion de los problemas que resuelve el sistema

Parallel                   Implementacion paralela que reusa todo el codigo secuencial
ParallelMethods            Implementacion de los distintos metodos de DVS en paralelo 
ParallelSchemeME           Implementacion del Esquema Maestro-Esclavo en paralelo

ParallelImplementation     Implementacion paralela usando el equivalente de ligas simbolicas al codigo
SequentialImplementation   Implementacion secuencial usando el equivalente de ligas simbolicas al codigo

ClassHierarchy             Contiene lo necesario para generar la jerarquia de clases en HTML, PDF y XML


Para compilar y correr la aplicacion segun sea el caso usar:

+ Codigo secuencial:
  $ cd SequentialImplementation
  $ ./Link
  $ make deps
  $ make 
  $ make run
o bien usar
  $ ./DVS-DDM file SimetricoExpXY -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
  $ make clean
  $ ./UnLink
  
+ Codigo paralelo:
  $ cd ParallelImplementation
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
  $ make clean
  $ ./UnLink
  
También es posible usar la compilación indicada en cada caso sin usar makefile 
especificada en README.txt del directorio correspondiente.


Parametros usados en la ejecución de DVS-DDM

Problema:
   file <archivo>


Malla:

   Malla gruesa:
      -Nx  <val = 3>
      -Ny  <val = Nx>
      -Nz  <val = Nx>

   Malla en cada subdominio
      -Mx  <val = Nx>
      -My  <val = Mx>
      -Mz  <val = Mx>


Método
   -m  <MF1 | PMF1 | MF2 | PMF2 | LM1 | PLM1 | LM2 | PLM2>
   

Visualizacion de resultados
   -p  <0 | 1 | 2 | 3>
        ^   ^   ^   ^
        |   |   |   |______visualiza todos los nodos x,y[,z],valor
        |   |   |__________visualiza por subdominio todos los nodos x,y[,z],valor
        |   |______________visualiza solo nodos de frontera interior x,y[,z],valor
        |__________________no visualiza

