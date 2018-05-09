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


#ifndef __LM1__
#define __LM1__

#include "Definiciones.hpp"
#include "DPMethod.hpp"
#include "CGM.hpp"
#include "DQGMRES.hpp"


class LM1 : public  DPMethod
{
private:
   ldouble *up;
   ldouble *uf;
   ldouble *f;
   EllipOp *op;


public:

   LM1(PropDef &props, EllipOp &op) : DPMethod(props)
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
