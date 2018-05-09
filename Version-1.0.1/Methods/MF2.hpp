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


#ifndef __MF2__
#define __MF2__


#include "Definiciones.hpp"
#include "DPMethod.hpp"
#include "PropDef.hpp"
#include "CGM.hpp"
#include "DQGMRES.hpp"



class MF2 : public  DPMethod
{
private:
   ldouble *up;
   EllipOp *op;

public:

   MF2(PropDef &props, EllipOp &op) : DPMethod(props)
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
