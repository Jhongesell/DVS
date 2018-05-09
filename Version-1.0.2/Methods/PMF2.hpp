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


#ifndef __PMF2__
#define __PMF2__


#include "Definiciones.hpp"
#include "DPMethod.hpp"
#include "CGM.hpp"
#include "DQGMRES.hpp"


class PMF2 : public  DPMethod
{
private:

   ldouble *up;
   ldouble *uf;
   EllipOp *op;

public:

   PMF2(PropDef &props, EllipOp &op) : DPMethod(props)
   {
      this->op = &op;
   }

   virtual void clean(void)
   {
      delete inter;
      inter = NULL;

      solver->clean();
      delete solver;
      solver = NULL;
   }

   ldouble dot(ldouble *u, ldouble *v);

   void multOp(ldouble *u, ldouble *v);

   void rhs(void);

   void solve(void);

};

#endif

