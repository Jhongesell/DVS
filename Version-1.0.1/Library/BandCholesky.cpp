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
#include <stdlib.h>
#include "Definiciones.hpp"
#include "BandCholesky.hpp"


void BandCholesky::convertBand(int n, MatrizDispersa *A)
{
   ce.nameClassFunct("BandCholesky", "convertBand");

   /* returns a banded Cholesky matrix from a square symmetrix matrix */
   int i, j, m, xcol, k, ind;
   bw = 0;
   for (i = 0; i < n; i++)
   {
      xcol = A->retornaNumeroColumnasBanda(i);
      for (ind = 0; ind < xcol; ind++)
      {
         k = A->retornaNumeroColumna(i, ind);
         m = (i >= k ? i - k : k - i);
         if (m > bw) bw = m;
      }
   }


//~ printf("\n\n banda %d \n",bw);
   AK = new (nothrow) ldouble *[bw + 1];
   if (AK == NULL) ce.memoryError("AK");
   for (i = 0; i < (bw + 1); i++)
   {
      AK[i] = new (nothrow) ldouble[n];
      if (AK[i] == NULL) ce.memoryError("AK[i]", i);
   }
   for (i = 0; i < (bw + 1); i++)
      for (j = 0; j < n; j++) AK[i][j] = 0.0;

   for (i = 0; i < n; i++)
   {
      xcol = A->retornaNumeroColumnasBanda(i);
      for (ind = 0; ind < xcol; ind++)
      {
         k = A->retornaNumeroColumna(i, ind);
         if (k <= i) AK[bw + k - i][i] = A->retornaValorColumna(i, ind);
      }
   }
}

void BandCholesky::factorLU(void)
{
   ldouble p;
   int i, j, k, m, m1;
   for (i = 0; i < n; i++)
   {
      m = (i + bw >= n ? n - 1 : i + bw);
      for (j = i + 1; j <= m; j++)
      {
         p = AK[bw + i - j][j] / AK[bw][i];
         m1 = (i + bw >= n ? n - 1 : i + bw);
         for (k = j; k <= m1; k++)
            AK[bw + j - k][k] -= p * AK[bw + i - k][k];
      }
      p = sqrt(AK[bw][i]);
      m = (i + bw >= n ? n - 1 : i + bw);
      for (k = i; k <= m; k++) AK[bw + i - k][k] /= p;
   }
}

void BandCholesky::solve(ldouble *x, ldouble *y)
{
   ldouble p;
   int i, j;
   for (i = 0; i < n; i++)
   {
      p = y[i];
      for (j = (i - bw < 0 ? 0 : i - bw); j < i; j++)
         p -= AK[bw + j - i][i] * x[j];
      x[i] = p / AK[bw][i];
   }
   for (i = n - 1; i >= 0; i--)
   {
      p = x[i];
      for (j = (i + bw >= n ? n - 1 : i + bw); j > i; j--)
         p -= AK[bw + i - j][j] * x[j];
      x[i] = p / AK[bw][i];
   }
}

void BandCholesky::print(void)
{
   int i, j, k;
   ldouble a;
   //~ printf("\nTamano %d, Banda %d\n",n, bw + 1);
   for (i = 0; i < n; i++)
   {
      printf("%d ", i);
      for (j = 0; j < bw + 1; j++)
      {
#ifdef __Double__
         printf(" %f ", AK[j][i] );
#else
         printf(" %Lf ", AK[j][i] );
#endif
      }
      printf("\n");
   }


   for (i = 0; i < n; i++)
   {
      printf("%d ", i);
      for (j = 0; j < n; j++)
      {
         k = i - j;
         a = (k < 0 ? (k < -bw ? 0.0 : AK[bw + k][j]) : (k > bw ? 0.0 : AK[bw - k][i]));
#ifdef __Double__
         printf(" %f ", a);
#else
         printf(" %Lf ", a);
#endif
      }
      printf("\n");
   }
}
