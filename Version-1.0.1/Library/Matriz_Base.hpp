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



#ifndef __Matriz_Base__
#define __Matriz_Base__

// Activar para depurar el codigo, desactivar para optimizar la velocidad de ejecucion
//#define DEPURAR



#include <string.h>
#include "ErrorControl.hpp"



/// Clase base para el trabajar con matrices
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
*/
class Matriz_Base
{

protected:
   /// Numero de columnas
   int Col;
   /// Numero de renglones
   int Ren;
   /// Tamano de la banda (solo si es bandada o dispersa)
   int Ban;
   /// Nombre de la matriz
   char *Nmb;
   /// Control de errores
   ErrorControl ce;

public:

   /// Constructor de la clase
   Matriz_Base(void)
   {
      Nmb = NULL;
      Col = Ren = Ban;
   }

   /// Destructor de la clase
   ~Matriz_Base()
   {
      delete []Nmb;
      Nmb = NULL;
   }

   /// Asigna nombre a la matriz
   /** @param nmb Nombre de la matriz */
   void asignaNombre(const char *nmb)
   {
      ce.nameClassFunct("Matriz_Base", "asignaNombre");

      Nmb = new (nothrow) char[strlen(nmb) + 2];
      if (Nmb == NULL) ce.memoryError("Nmb");
      strcpy(Nmb, nmb);
   }

   /// Retorna el numero de renglones de la matriz
   /** @return Regresa el numero de renglones de la matriz */
   inline int renglones(void)
   {
      return Ren;
   }

   /// Retorna el numero de columnas de la matriz
   /** @return Regresa el numero de columnas de la matriz */
   inline int columnas(void)
   {
      return Col;
   }


};

#endif
