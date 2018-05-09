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


#ifndef __Solvable__
#define __Solvable__

#include "Definiciones.hpp"
#include <stdlib.h>
#include <stdio.h>


class Solvable
{

protected:

   const char *name;

public:

   Solvable(void)
   {
      name = NULL;
   }

   virtual ~Solvable(void)
   {
   }

   virtual void clean(void) = 0;

   virtual void solve(ldouble *x, ldouble *y) = 0;

   virtual int getIter(void) = 0;

   const char *getName(void)
   {
      return name;
   }

};

#endif
