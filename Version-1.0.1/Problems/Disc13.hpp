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


#ifndef  __Disc13__
#define __Disc13__


#include <math.h>
#include "Definiciones.hpp"
#include "FunctionV1.hpp"



class Disc13 : public FunctionV1
{
private:

   ldouble var;


public:

   Disc13(ldouble b)
   {
      var = b;
   }

   // Condiciones de Frontera Tadeus
   inline ldouble eval(int d, ldouble *x)
   {

      if (fabs(x[0] - -1.0) < EPS_EQUAL && x[1] > -1.0 && x[1] <= 1.0)  return 1.0;
      if (fabs(x[1] - 1.0) < EPS_EQUAL && x[0] >= -1.0 && x[0] <= 1.0)  return 1.0;

      if (fabs(x[1] - -1.0) < EPS_EQUAL && x[0] >= -1.0 && x[0] <= 1.0) return 0.0;

      if (fabs(x[0] - 1.0) < EPS_EQUAL  && x[1] >= -1.0 && x[1] <= 1.0)  return ( (1.0 + x[1]) / 2.0);


      return 0.0;
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
