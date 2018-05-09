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


#ifndef __FunctionV__
#define __FunctionV__

#include "Definiciones.hpp"


class FunctionV
{
protected:
   int dim;


public:

   FunctionV(void)
   {
      dim = 0;
   }

   virtual ~FunctionV()
   {}

   virtual ldouble eval(int d, ldouble *x) = 0;

   void dimension(int d)
   {
      dim = d;
   }

};

#endif
