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


#ifndef __ExpVXYZ__
#define __ExpVXYZ__


#include <math.h>
#include "Definiciones.hpp"
#include "FunctionV1.hpp"


class ExpVXYZ : public FunctionV1
{
private:

   ldouble var;


public:

   ExpVXYZ(ldouble b)
   {
      var = b;
   }

   inline ldouble eval(int d, ldouble *x)
   {
      int i;
      ldouble v[] = {1.0, 1.0, 1.0};
      ldouble val = 0.0;
      for (i = 0; i < 3; i++) val += x[i] * v[i];
      return var * exp(val);
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
