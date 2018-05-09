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



#ifndef __AllPrimal__
#define __AllPrimal__

#include "Primal.hpp"
#include <string>


class AllPrimal : public Primal
{

public:

   const char *name;

   AllPrimal(void) : name("AllPrimal")
   {  }

   inline bool isPrimal(int type, int *coordN, int *coordM)
   {
      //System.out.println(name);
      return (type & INTBD) != 0;
   }

};

#endif


