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


#ifndef __ExpX__
#define __ExpX__


#include <math.h>
#include "Definiciones.hpp"
#include "FunctionV1.hpp"


class ExpX : public FunctionV1
{
private:

   ldouble var;


public:

   ExpX(ldouble b)
   {
      var = b;
   }

   inline ldouble eval(int d, ldouble *x)
   {
      return var * exp(x[0]);
   }

   inline ldouble getVar(void)
   {
      return var;
   }

   inline void setVar(ldouble b)
   {
      var = b;
   }

};

#endif
