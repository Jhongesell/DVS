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

Si se quiere usar make
  $ ./Link
  $ make deps
  $ make 
  $ make run
  $ make clean
  $ ./UnLink 
o bien usar, por ejemplo
  $ ./Link
  $ make deps
  $ make 
  $ ./DVS-DDM file prob1 -m MF1  -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
  $ make clean
  $ ./UnLink


Si no se quiere usar make, para activar los fuentes respectivos
  $ ./Link
Compilar usando:
  $ g++ -O1 *.cpp
Ejecutar en Linux:
  $ ./a.out file prob1 -m MF1 -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10

Para hacer este analisis de rendimiento, hacer:
  $ g++ -g -pg *.cpp
  $ ./a.out file prob1 -m MF1 -p 3 -Nx 10 -Ny 4 -Mx 10 -My 10
  $ gprof -c -z a.out > sal.txt
El archivo sal.txt contiene el analisis de rendimiento detallado.

Para rastreo de problemas con la manipulacion de memoria y punteros desbordados:
  $ g++ -g *.cpp
  $ valgrind --tool=memcheck --leak-check=yes --show-reachable=yes ./a.out file prob7 -m MF1 -Nx 3
