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


#ifndef __Constant__
#define __Constant__

#include "Definiciones.hpp"
#include "FunctionV1.hpp"


class Constant : public FunctionV1
{
private:

   ldouble a;


public:

   Constant(ldouble b)
   {
      a = b;
   }

   inline ldouble eval(int d, ldouble *x)
   {
      return a;
   }

   inline ldouble getVar(void)
   {
      return a;
   }

   inline void setVar(ldouble b)
   {
      a = b;
   }

};

#endif
