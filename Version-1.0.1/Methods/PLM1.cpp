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
#include "PLM1.hpp"


// DVS-FETI-DP

ldouble PLM1::dot(ldouble *u, ldouble *v)
{
   ldouble val = 0.0;
   dualp->solveS(scr, u);
   for (int i = 0; i < nDual; i++) val += scr[i] * v[i];
   return val;
}

void PLM1::multOp(ldouble *u, ldouble *v)
{
   dualp->solveS(uf, u);
   dualp->j(uf, v);
   dualp->multS(v, uf);
   dualp->j(uf, v);
}

void PLM1::rhs(void)
{
   inter->rhs(0);
   dualp->solveAPP(0, 0, 1, 2);
   inter->multOp(1, 1);
   inter->diff(0, 0, 1);
   dualp->fromSubdomains(0, f);

   dualp->solveS(up, f);
   dualp->j(up, rhss);

   dualp->multS(rhss, uf);
   dualp->j(uf, rhss);
}


void PLM1::solve(void)
{
   ce.nameClassFunct("PLM1", "solve");
   DPMethod::initialize();

   int i;
   scr = new (nothrow) ldouble[nDual];
   if (scr == NULL) ce.memoryError("scr");
   up = new (nothrow) ldouble[nDual];
   if (up == NULL) ce.memoryError("up");
   uf = new (nothrow) ldouble[nDual];
   if (uf == NULL) ce.memoryError("uf");
   f = new (nothrow) ldouble[nDual];
   if (f == NULL) ce.memoryError("f");
   for (i = 0; i < nDual; i++) scr[i] = up[i] = uf[i] = f[i] = 0.0;


   MultOp *pt = (MultOp*) this;
   DotProd *ptd = (DotProd*) this;
   int n = (int) sqrt(pt->getSize());

   if (op->isSymmetric()) solver = new (nothrow) CGM(*pt, *ptd, epsilon);
   else solver = new (nothrow) DQGMRES(*pt, n, epsilon);
   if (solver == NULL) ce.memoryError("solver");

   rhs();

   solver->solve(u, rhss);

   for (i = 0; i < nDual; i++) u[i] = f[i] - u[i];
   dualp->solveS(up, u);
   dualp->a(up, uf);

   dualp->calcValues(uf);
   print(uf);

#ifdef NUMERO_CONDICIONAMIENTO
   conditionalNumber(op->isSymmetric());
#endif

   delete []scr;
   delete []up;
   delete []uf;
   delete []f;
}
