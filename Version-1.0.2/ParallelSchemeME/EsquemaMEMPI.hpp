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


#ifndef __EsquemaMEMPI__
#define __EsquemaMEMPI__


//#define SUBDOMINIOS_CONTIGUOS


#include "mpi.h"
#include "ErrorControl.hpp"



/// Clase base para definir el Esquema Maestro-Esclavo en MPI
/** Clase base para definir el Esquema Maestro-Esclavo para
    programar en paralelo mediante el paso de mensajes usando
    MPI, donde el primer procesador (id = 0) es el nodo mestro
    y el resto son los nodos esclavos. Las tareas se pueden repartir
    de manara que subdominios contiguos queden en un mismo
    nodo esclavo o queden en distinto nodo esclavo.
*/
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.0
    @bug No hay errores conocidos
*/
class EsquemaMEMPI
{

protected:

   /// Identificador
   int id;
   /// Numero de procesadores
   int np;
   /// Numero de tareas por nodo esclavo
   int *ta;
   /// Numero de nodos esclavos a utilizar (los que tienen carga)
   int npu;
   /// Control de errores
   ErrorControl ce;


public:

   /// Constructor de la clase
   /** @param id Identificador
       @param np Numero de procesadores */
   EsquemaMEMPI(int id, int np)
   {
      this->id = id;
      this->np = np;
      ta = NULL;
   }

   /// Destructor de la clase
   ~EsquemaMEMPI()
   {
      if (ta) delete []ta;
   }

   /// Genera el reparto de carga
   /** @param n Numero de trabajos */
   void generaRepartoCarga(int n);

   /// Reparte la carga de trabajo entre los nodos esclavos
   /** @param np Numero de procesador esclavo
       @param st  Indice de tarea dentro del nodo esclavo
       @param tarea Tarea la cual debe ser repartida   */
   void reparteCargaTrabajo(int &np, int &ind, int tarea);

   /// Retorna el numero de procesadores a usar por el esquema M-E
   /** @return  Numero de procesadores a usar dentro del esquema Maestro-Esclavo */
   inline int numeroProcesadoresUsar(void)
   {
      return npu;
   }

};

#endif

/*
Reparto de carga en el esquema Maestro-Esclavo
+ Asignacion de subdominios contiguos en un mismo nodo esclavo
   - Permite hacer una asignacion o recoleccion de datos en bloque
+ Asignacion de subdominios alternados entre nodos esclavos


Por considerar
+Poder configurar cuando un nodo es mas lento que los otros para hacer
   heterogenea la reparticion de carga

*/
