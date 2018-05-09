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


#ifndef __SinPixSinPiy__
#define __SinPixSinPiy__


#include <math.h>
#include "Definiciones.hpp"
#include "FunctionV1.hpp"



class SinPixSinPiy : public FunctionV1
{
private:

   ldouble var;


public:

   SinPixSinPiy(ldouble b)
   {
      var = b;
   }

   inline ldouble eval(int d, ldouble *x)
   {
      return var * sin(M_PI * x[0]) * sin(M_PI * x[1]);
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
