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
#include "MultBandSym.hpp"


void MultBandSym::multOp(ldouble *x, ldouble *y)
{
   int i, j, m;
   for (i = 0; i < n; i++)
   {
      y[i] = 0.0;
      for (j = (i - bw < 0 ? 0 : i - bw); j <= i; j++)
         y[i] += AK[bw + j - i][i] * x[j];
      m = (i + bw >= n ? n - 1 : i + bw);
      for (j = i + 1; j <= m; j++)
         y[i] += AK[bw + i - j][j] * x[j];
   }
}
