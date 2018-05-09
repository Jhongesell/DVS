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



#ifndef __EllipOp__
#define __EllipOp__

#include "Definiciones.hpp"
#include "FunctionV.hpp"
#include <stdlib.h>


/* Represents the constant coefficient elliptic operator
   -a(del^2)u + (b dot grad(u)) + cu = f(x) with Dirichlet boundary
   condition u(x) = g(x) on the external boundary
*/
class EllipOp
{

public:

   int nDim;
   ldouble *a;
   ldouble *b;
   ldouble c;
   FunctionV *f;     // r.h.s. of operator
   FunctionV *g;     // Dirichlet boundary condition
   FunctionV *sol;   // Problem solution id known


   EllipOp(int nDim, ldouble *a, ldouble *b, ldouble c, FunctionV &f, FunctionV &g, FunctionV &sol)
   {
      this->nDim = nDim;
      //~ if (b.size() == 0) b.resize(3);
      this->a = a;
      this->b = b;
      this->c = c;
      this->f = &f;
      this->g = &g;
      this->sol = &sol;
   }


   EllipOp(int nDim, ldouble *a, ldouble *b, ldouble c)
   {
      this->nDim = nDim;
      this->a = a;
      this->b = b;
      this->c = c;
      this->f = NULL;
      this->g = NULL;
      this->sol = NULL;
   }

   ~EllipOp()
   {
      f = NULL;
      g = NULL;
      sol = NULL;
      a = NULL;
      b = NULL;
   }

   ldouble *getA(void)
   {
      return a;
   }

   ldouble *getB(void)
   {
      return b;
   }

   ldouble getC(void)
   {
      return c;
   }

   FunctionV *getF(void)
   {
      return f;
   }

   FunctionV *getG(void)
   {
      return g;
   }

   void setF(FunctionV &f)
   {
      this->f = &f;
   }

   void setG(FunctionV &g)
   {
      this->g = &g;
   }

   bool isSymmetric(void)
   {
      if (b == NULL || (b[0] == 0.0 && b[1] == 0.0 && b[2] == 0.0)) return true;
      return false;
   }

};

#endif
