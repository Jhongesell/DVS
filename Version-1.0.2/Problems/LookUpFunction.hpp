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



#ifndef __LookUpFunction__
#define __LookUpFunction__


#include <string>
#include "FunctionV1.hpp"
#include "ErrorControl.hpp"



class LookUpFunction
{
protected:
   /// Control de errores
   ErrorControl ce;


public:

   FunctionV1 *getF(char *s);


};

#endif

