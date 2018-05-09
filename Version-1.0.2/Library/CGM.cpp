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


#include <math.h>
#include "CGM.hpp"


void CGM::solve(ldouble *u, ldouble *b)
{
   int i;
   nIter = 1;
   A->multOp(u, v);    // v = A*u
   for (i = 0; i < n; i++) p[i] = r[i] = b[i] - v[i];
   gamma = dotP->dot(r, r);
   mu = norm(r);
   while (mu > eps && nIter < nMaxIter)
   {
#ifdef RESIDUAL
      printf("\nIterations:  %d, Error: %1.16f", nIter, mu);
      fflush(stdout);
#endif
      A->multOp(p, v);    //  v = A*p
      lambda = dotP->dot(p, v);   // lambda = p*A*p
      alpha = gamma / lambda; // alpha = r*r/p*A*p
      for (i = 0; i < n; i++)
      {
         u[i] += alpha * p[i];
         r[i] -= alpha * v[i];
      }
      lambda = dotP->dot(r, r);
      beta = lambda / gamma;
      for (i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
      gamma = lambda;
      mu = norm(r);
      nIter++;
   }
#ifdef RESIDUAL
   printf("\nIterations:  %d, Error: %1.16f", nIter, mu);
   fflush(stdout);
#endif
}

ldouble CGM::norm(ldouble *x)
{
   ldouble v, val = 0.0;
   for (int i = 0; i < n; i++)
      if ((v = fabs(x[i])) > val) val = v;
   return val;
}

