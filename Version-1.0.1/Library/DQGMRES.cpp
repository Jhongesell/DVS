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
#include "DQGMRES.hpp"



void DQGMRES::inicializa(void)
{
   ce.nameClassFunct("DQGMRES", "inicializa");

   int i, j;
   k1 = k + 1;

   p = new (nothrow) ldouble *[k1];
   if (p == NULL) ce.memoryError("p");
   for (i = 0; i < k1; i++)
   {
      p[i] = new (nothrow) ldouble[n];
      if (p[i] == NULL) ce.memoryError("p[i]", i);
   }
   for (i = 0; i < k1; i++)
      for (j = 0; j < n; j++) p[i][j] = 0.0;


   cs = new (nothrow) ldouble *[k1];
   if (cs == NULL) ce.memoryError("cs");
   for (i = 0; i < k1; i++)
   {
      cs[i] = new (nothrow) ldouble[2];
      if (cs[i] == NULL) ce.memoryError("cs[i]", i);
   }
   for (i = 0; i < k1; i++)
      for (j = 0; j < 2; j++) cs[i][j] = 0.0;

   h = new (nothrow) ldouble *[k1 + 1];
   if (h == NULL) ce.memoryError("h");
   for (i = 0; i < (k1 + 1); i++)
   {
      h[i] = new (nothrow) ldouble[k1];
      if (h[i] == NULL) ce.memoryError("h[i]", i);
   }
   for (i = 0; i < (k1 + 1); i++)
      for (j = 0; j < k1; j++) h[i][j] = 0.0;

   q = new (nothrow) ldouble *[k1 + 1];
   if (q == NULL) ce.memoryError("q");
   for (i = 0; i < (k1 + 1); i++)
   {
      q[i] = new (nothrow) ldouble[n];
      if (q[i] == NULL) ce.memoryError("q[i]", i);
   }
   for (i = 0; i < (k1 + 1); i++)
      for (j = 0; j < n; j++) q[i][j] = 0.0;

   v = new (nothrow) ldouble [n];
   if (v == NULL) ce.memoryError("v");
   for (i = 0; i < n; i++) v[i] = 0.0;
}


void DQGMRES::applyOmega(int m)
{
   int i, i1, m1 = m % k1;
   ldouble z1, z2;
   for (i = (m - k < 0 ? 0 : m - k); i < m; i++)
   {
      i1 = i % k1;
      z1 = cs[i1][0] * h[i1][m1] + cs[i1][1] * h[i1 + 1][m1];
      z2 = -cs[i1][1] * h[i1][m1] + cs[i1][0] * h[i1 + 1][m1];
      h[i1][m1] = z1;
      h[i1 + 1][m1] = z2;
   }
}



void DQGMRES::solve(ldouble *x, ldouble *b)
{
   int i, i1, j, j1, m, m1 = 0;
   ldouble val = 0;
   gm = 0.0;
   nIter = 1;
   for (i = 0; i < n; i++)
   {
      x[i] = 0.0;
      gm += b[i] * b[i];
   }

   if (gm == 0.0) return;
   gm = sqrt(gm);
   for (i = 0; i < n; i++) q[0][i] = b[i] / gm;
   for (m = 0; m < nMaxIter; m++)
   {
      m1 = m % k1;
      mult->multOp(q[m1], v);
      for (i = (m - k + 1 < 0 ? 0 : m - k + 1); i <= m; i++)
      {
         i1 = i % k1;
         val = 0.0;
         for (j = 0; j < n; j++) val += v[j] * q[i1][j];
         h[i1][m1] = val;
         for (j = 0; j < n; j++) v[j] = v[j] - val * q[i1][j];
      }
      val = 0.0;
      for (j = 0; j < n; j++) val += v[j] * v[j];
      val = sqrt(val);
      h[m1 + 1][m1] = val;
      for (j = 0; j < n; j++) q[m1 + 1][j]  = v[j] / val;
      applyOmega(m);
      val = sqrt(h[m1][m1] * h[m1][m1] + h[m1 + 1][m1] * h[m1 + 1][m1]);
      cs[m1][0] = h[m1][m1] / val;
      cs[m1][1] = h[m1 + 1][m1] / val;
      gm1 = -cs[m1][1] * gm;
      gm = cs[m1][0] * gm;
      h[m1][m1] = cs[m1][0] * h[m1][m1] + cs[m1][1] * h[m1 + 1][m1];
      for (i = 0; i < n; i++)
      {
         val = q[m1][i];
         for (j = (m - k < 0 ? 0 : m - k); j < m; j++)
         {
            j1 = j % k1;
            val -= h[j1][m1] * p[j1][i];
         }
         p[m1][i] = val / h[m1][m1];
         x[i] += gm * p[m1][i];
      }
      val = (gm1 < 0.0 ? -gm1 : gm1);
      gm = gm1;
#ifdef RESIDUAL
      printf("\nIterations:  %d, Error: %1.16f gm %1.16f ", nIter, val, gm);
      fflush(stdout);
#endif
      if (val < eps) break;
      nIter++;
   }
}

