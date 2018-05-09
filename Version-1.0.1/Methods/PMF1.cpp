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
#include "PMF1.hpp"




ldouble PMF1::dot(ldouble *u, ldouble *v)
{
   ldouble val = 0.0;
   dualp->multS(u, scr);
   for (int i = 0; i < nDual; i++) val += scr[i] * v[i];
   return val;
}

void PMF1::multOp(ldouble *u, ldouble *v)
{
   dualp->multS(u, v);
   dualp->a(v, v);
   dualp->solveS(v, v);
   dualp->a(v, v);
}

void PMF1::rhs(void)
{

   inter->rhs(0);
   dualp->solveAPP(0, 0, 1, 2);
   inter->multOp(1, 1);
   inter->diff(0, 0, 1);

   dualp->fromSubdomains(0, rhss);
   dualp->a(rhss, rhss);

   dualp->solveS(rhss, rhss);
   dualp->a(rhss, rhss);
}

void PMF1::solve(void)
{
   ce.nameClassFunct("PMF1", "solve");
   DPMethod::initialize();

   int i;
   scr = new (nothrow) ldouble[nDual];
   if (scr == NULL) ce.memoryError("scr");
   for (i = 0; i < nDual; i++) scr[i] = 0.0;

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

#ifdef NUMERO_CONDICIONAMIENTO
   conditionalNumber(op->isSymmetric());
#endif

   delete []scr;
}
