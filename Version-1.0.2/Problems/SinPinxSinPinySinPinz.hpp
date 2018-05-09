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


#ifndef __SinPinxSinPinySinPinz__
#define __SinPinxSinPinySinPinz__


#include <math.h>
#include "Definiciones.hpp"
#include "FunctionV1.hpp"



class SinPinxSinPinySinPinz : public FunctionV1
{
private:

   ldouble var;
   ldouble n;

public:

   SinPinxSinPinySinPinz(ldouble b)
   {
      var = b;
      n = 4.0;
   }

   inline ldouble eval(int d, ldouble *x)
   {
      return 3.0 * n * n * M_PI * M_PI * sin(M_PI * n * x[0]) * sin(M_PI * n * x[1]) * sin(M_PI * n * x[2]);
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
