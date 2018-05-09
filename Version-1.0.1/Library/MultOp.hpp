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



#ifndef __MultOp__
#define __MultOp__

#include "Definiciones.hpp"


class MultOp
{

public:

   /// y = A*x
   virtual void multOp(ldouble *x, ldouble *y) = 0;

   /// vector size
   virtual int getSize(void) = 0;

};

#endif
