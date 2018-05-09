// Grupo de Modelacion Matematica y Computacional
// Instituto de Goefisica
// Universidad Nacional Autonoma de Mexico
// Mexico D.F.
//
// Autores:
// Ismael Herrera Revilla
// Robert Yates Smit
// Antonio Carrillo Ledesma
//
// Codigo liberado bajo la licencia GPL ver. 2.0


#include "EsquemaMEMPI.hpp"

/*
Reparto de carga en el esquema Maestro-Esclavo
+ Asignacion de subdominios contiguos en un mismo nodo esclavo
   - Permite hacer una asignacion o recoleccion de datos en bloque
+ Asignacion de subdominios alternados entre nodos esclavos


Por considerar
+Poder configurar cuando un nodo es mas lento que los otros para hacer
   heterogenea la reparticion de carga

*/



// Genera el reparto de carga
// @param n Numero de trabajos
void EsquemaMEMPI::generaRepartoCarga(int n)
{
   ce.nameClassFunct("EsquemaMEMPI", "generaRepartoCarga");
   int i, j;

   // Ajusta el numero de procesadores a usar
   if (n >= (np - 1)) npu = np - 1;
   else npu = n;

   // Arreglo para soportar el reparto de carga por cada esclavo
   ta = new (nothrow) int [npu];
   if (ta == NULL) ce.memoryError("ta");
   for (i = 0; i < npu; i++) ta[i] = 0;

   // Reparto de carga
   j = 0, i = 0;
   while (i < n)
   {
      ta[j] ++;
      i++;
      j++;
      if (j >= npu) j = 0;
   }

   // Visualiza la carga
   for (i = 0; i < npu; i++) printf("\nId = %d  ==> TAREAS = %d ", i + 1, ta[i]);
   printf("\n\nNumero de procesadores esclavos a usar = %d\n\n", npu);
   fflush(stdout);

   // Avisa a cada esclavo la cantidad de tareas a soportar
   for (i = 1; i <= npu; i++)  MPI::COMM_WORLD.Send(&ta[i - 1], 1, MPI::INT, i, 1);
   j = 0;
   for (i = (npu + 1); i < np; i++) MPI::COMM_WORLD.Send(&j, 1, MPI::INT, i, 1);
}


#ifdef SUBDOMINIOS_CONTIGUOS
// El reparto de carga se hace de forma que nodos contiguos queden en un mismo nodo esclavo
/// Reparte la carga de trabajo entre los nodos esclavos
/** @param np Numero de procesador esclavo
    @param st  Indice de tarea dentro del nodo esclavo
    @param tarea Tarea la cual debe ser repartida   */
void EsquemaMEMPI::reparteCargaTrabajo(int &np, int &ind, int tarea)
{
   int i = 0;
   // Indica el numero de tarea dentro del esclavo
   ind = tarea;
   while (ind >= ta[i])
   {
      ind -= ta[i];
      i++;
   }
   // Indica el numero de esclavo en el que estara la tarea
   np = i + 1;
}

#else

/// Reparte la carga de trabajo entre los nodos esclavos
// El reparto de carga entre los nodos esclavos, si hay mas tareas que procesadores se asignaran mas de una tarea a cada esclavo, el reparto no necesariamente es homogeneo
/** @param np Numero de procesador esclavo
    @param st  Indice de tarea dentro del nodo esclavo
    @param tarea Tarea la cual debe ser repartida   */
void EsquemaMEMPI::reparteCargaTrabajo(int &np, int &ind, int tarea)
{
   // Indica el numero de esclavo en el que estara la tarea
   np = (tarea % npu) + 1;
   // Indica el numero de tarea dentro del esclavo
   ind = tarea / npu;
}
#endif


