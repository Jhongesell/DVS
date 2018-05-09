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



#ifndef __NoPrimal__
#define __NoPrimal__

#include "Primal.hpp"


class NoPrimal : public Primal
{

public:

   const char *name;

   NoPrimal(void) : name("NoPrimal")
   {  }


   bool isPrimal(int type, int *coordN, int *coordM)
   {
      //printf("%s",name);
      return false;
   }

};

#endif
