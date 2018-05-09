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



#ifndef __FunctionV1__
#define __FunctionV1__

#include "Definiciones.hpp"
#include "FunctionV.hpp"


class FunctionV1 : public FunctionV
{

public:

   virtual void setVar(ldouble x) = 0;
   virtual ldouble getVar(void) = 0;

};

#endif
