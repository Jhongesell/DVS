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


#include "Definiciones.hpp"
#include "MF1.hpp"


ldouble MF1::dot(ldouble *u, ldouble *v)
{
   ldouble val = 0.0;
   for (int i = 0; i < nDual; i++) val += u[i] * v[i];
   return val;
}

void MF1::multOp(ldouble *u, ldouble *v)
{
   dualp->multS(u, v);
   dualp->a(v, v);
}

void MF1::rhs(void)
{
   inter->rhs(0);
   dualp->solveAPP(0, 0, 1, 2);
   inter->multOp(1, 1);
   inter->diff(0, 0, 1);

   dualp->fromSubdomains(0, rhss);
   dualp->a(rhss, rhss);
}

void MF1::solve(void)
{
   ce.nameClassFunct("MF1", "solve");
   DPMethod::initialize();

   MultOp *pt = (MultOp*) this;
   DotProd *ptd = (DotProd*) this;
   int n = (int) sqrt(pt->getSize());

   if (op->isSymmetric()) solver = new (nothrow) CGM(*pt, *ptd, epsilon);
   else solver = new (nothrow) DQGMRES(*pt, n, epsilon);
   if (solver == NULL) ce.memoryError("solver");

   rhs();
   solver->solve(u, rhss);
   dualp->calcValues(u);
   print(u);
}
