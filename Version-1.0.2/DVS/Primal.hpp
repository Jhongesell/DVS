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


#ifndef __Primal__
#define __Primal__


class Primal
{

public:

   Primal()
   {}

   virtual ~Primal()
   {}

   static const int KNOWN = 1, INTERIOR = 2, INTBD = 4, VERTEX = 8, EDGE = 16, FACE = 32, PRIMAL = 64;

   virtual bool isPrimal(int type, int *coordN, int *coordM) = 0;

};

#endif
