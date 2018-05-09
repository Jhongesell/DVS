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
#include "PMF2.hpp"




ldouble PMF2::dot(ldouble *u, ldouble *v)
{
   ldouble val = 0.0;
   dualp->multS(u, scr);
   for (int i = 0; i < nDual; i++) val += scr[i] * v[i];
   return val;
}

void PMF2::multOp(ldouble *u, ldouble *v)
{
   dualp->j(u, v);
   dualp->multS(v, uf);
   dualp->j(uf, v);
   dualp->solveS(v, v);
}

void PMF2::rhs(void)
{
   inter->rhs(0);
   dualp->solveAPP(0, 0, 1, 2);
   inter->multOp(1, 1);
   inter->diff(0, 0, 1);

   dualp->fromSubdomains(0, rhss);

   dualp->solveS(up, rhss);
   dualp->j(up, uf);

   dualp->multS(uf, rhss);
   dualp->j(rhss, rhss);
   dualp->solveS(rhss, rhss);
}


void PMF2::solve(void)
{
   ce.nameClassFunct("PMF2", "solve");
   DPMethod::initialize();

   int i;
   scr = new (nothrow) ldouble[nDual];
   if (scr == NULL) ce.memoryError("scr");
   up = new (nothrow) ldouble[nDual];
   if (up == NULL) ce.memoryError("up");
   uf = new (nothrow) ldouble[nDual];
   if (uf == NULL) ce.memoryError("uf");
   for (i = 0; i < nDual; i++) scr[i] = up[i] = uf[i] = 0.0;

   MultOp *pt = (MultOp*) this;
   DotProd *ptd = (DotProd*) this;
   int n = (int) sqrt(pt->getSize());

   if (op->isSymmetric()) solver = new (nothrow) CGM(*pt, *ptd, epsilon);
   else solver = new (nothrow) DQGMRES(*pt, n, epsilon);
   if (solver == NULL) ce.memoryError("solver");

   rhs();
   solver->solve(u, rhss);

   for (i = 0; i < nDual; i++) u[i] = up[i] - u[i];

   dualp->calcValues(u);
   print(u);

#ifdef NUMERO_CONDICIONAMIENTO
   conditionalNumber(op->isSymmetric());
#endif

   delete []scr;
   delete []up;
   delete []uf;
}

