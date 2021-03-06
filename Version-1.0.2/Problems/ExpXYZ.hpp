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


#ifndef  __ExpXYZ__
#define __ExpXYZ__


#include <math.h>
#include "Definiciones.hpp"
#include "FunctionV1.hpp"



class ExpXYZ : public FunctionV1
{
private:

   ldouble var;


public:

   ExpXYZ(ldouble b)
   {
      var = b;
   }

   inline ldouble eval(int d, ldouble *x)
   {
      return  exp(x[0] * x[1] * x[2]);
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
