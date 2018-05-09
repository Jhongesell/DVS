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



#ifndef __DotProd__
#define __DotProd__

#include "Definiciones.hpp"


class DotProd
{

public:

   // returns  x (dot) y
   virtual ldouble dot(ldouble *x, ldouble *y) = 0;

};

#endif
