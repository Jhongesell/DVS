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


#ifndef __NSfExpXY__
#define __NSfExpXY__

#include "FunctionV1.hpp"
#include "Definiciones.hpp"
#include <math.h>



class NSfExpXY : public FunctionV1
{
private:

   double var;


public:

   NSfExpXY(double b)
   {
      var = b;
   }

   inline double eval(int d, double *x)
   {
      double p = (var == 0.0 ? x[0] + x[1] - x[0] * x[0] - x[1] * x[1] : var); // Si var es igual a cero entonces p es igual a x[0]+ x[1]- x[0]*x[0] - x[1]*x[1]
      return p * exp(x[0] * x[1]);
   }

   inline double getVar(void)
   {
      return var;
   }

   inline void setVar(double b)
   {
      var = b;
   }

};

#endif
