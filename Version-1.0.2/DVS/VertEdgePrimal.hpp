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



#ifndef __VertEdgePrimal__
#define __VertEdgePrimal__

#include "Primal.hpp"


class VertEdgePrimal : public Primal
{

public:

   const char *name;

   VertEdgePrimal(void) : name("VertEdgePrimal")
   { }

   inline bool isPrimal(int type, int *coordN, int *coordM)
   {
      //printf("%s",name);
      return ((type & VERTEX) != 0 || (type & EDGE) != 0);
   }

};

#endif
