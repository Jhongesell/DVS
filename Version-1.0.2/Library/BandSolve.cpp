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



#include <stdlib.h>
#include "Definiciones.hpp"
#include "BandSolve.hpp"


/* Convert square matrix A into banded matrix and then do a standard BandSolve
   Leaves bandwidth in global variable bw  */
BandSolve::BandSolve(int n, MatrizDispersa *A)
{
   //printf("\nBandSolve %d",n);
   convertBand(n, A);
   this->n = n;
   factorLU();
}


/* Convert square matrix A into banded matrix and then do a standard BandSolve
   Leaves bandwidth in global variable bw  */
BandSolve::BandSolve(int n, ldouble **A)
{
   //printf("\nBandSolve %d",n);
   convertBand(n, A);
   this->n = n;
   factorLU();
}


void BandSolve::factorLU(void)
{
   /* Decomposes the banded matrix a as a = LU  -> stores the result over a */
   ldouble p, q;
   int i, j, k;
   for (i = 0; i < n; i++)
   {
      p = AK[bw][i];    // ith diagonal -> pivot
      for (j = i + 1; j < n && j <= bw + i; j++)
      {
         q = AK[bw + i - j][j] / p;
         AK[bw + i - j][j] = q;
         if (q != 0.)
            for (k = i + 1; k < n && k <= bw + i; k++)
               AK[bw + k - j][j] -= AK[bw + k - i][i] * q;
      }
   }
}

void BandSolve::solve(ldouble *x, ldouble *y)
{
   /* Solves the equation AK*x=y for x where a has been previously factored */
   ldouble p;
   int i, j;
   x[0] = y[0];
   for (i = 1; i < n; i++)
   {
      p = y[i];
      for (j = i - 1; j >= 0 && j >= i - bw; j--)
         p -= AK[bw + j - i][i] * x[j];
      x[i] = p;
   }
   x[n - 1] /= AK[bw][n - 1];
   for (i = n - 2; i >= 0; i--)
   {
      p = x[i];
      for (j = i + 1; j < n && j <= bw + i; j++)
         p -= AK[bw + j - i][i] * x[j];
      x[i] = p / AK[bw][i];
   }
}

void BandSolve::convertBand(int n, ldouble **A)
{
   ce.nameClassFunct("BandSolve", "convertBand");
   int i, j, m = 0;
   bw = 0;
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
      {
         m = (i >= j ? i - j : j - i);
         if (m > bw && A[i][j] != 0.0) bw = m;
      }

   int tm = 2 * bw + 1;
   //printf("\nbandwidth %d %d",n,(2*bw + 1));
   AK = new (nothrow) ldouble *[tm];
   if (AK == NULL) ce.memoryError("AK");
   for (i = 0; i < tm; i++)
   {
      AK[i] = new (nothrow) ldouble[n];
      if (AK[i] == NULL) ce.memoryError("AK[i]", i);
      for (j = 0; j < n; j++) AK[i][j] = 0.0;
   }

   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
      {
         if (A[i][j] != 0.0) AK[bw + j - i][i] = A[i][j];
      }
}

void BandSolve::convertBand(int n, MatrizDispersa *A)
{
   ce.nameClassFunct("BandSolve", "convertBand");

   int i, j, m = 0, xcol, k , ind;
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

   //printf("\nbandwidth %d %d",n,(2*bw + 1));
   AK = new (nothrow) ldouble *[2 * bw + 1];
   if (AK == NULL) ce.memoryError("AK");
   for (i = 0; i < (2 * bw + 1); i++)
   {
      AK[i] = new (nothrow) ldouble[n];
      if (AK[i] == NULL) ce.memoryError("AK[i]", i);
   }
   for (i = 0; i < (2 * bw + 1); i++)
      for (j = 0; j < n; j++) AK[i][j] = 0.0;

   for (i = 0; i < n; i++)
   {
      xcol = A->retornaNumeroColumnasBanda(i);
      for (ind = 0; ind < xcol; ind++)
      {
         k = A->retornaNumeroColumna(i, ind);
         AK[bw + k - i][i] = A->retornaValorColumna(i, ind);
      }
   }
}


void BandSolve::print(void)
{
   /* Prints the banded matrix to standard output */
   printf("\nBandSolve");
   int i, j;
   for (i = 0; i < n; i++)
   {
      printf("row  %d ", i);
      for (j = 0; j < n; j++)
         if (abs(i - j) > bw) printf(" 0.0 ");
#ifdef __Double__
         else printf(" %f ", AK[bw + j - i][i]);
#else
         else printf(" %Lf ", AK[bw + j - i][i]);
#endif
      printf("\n");
   }
   printf("\n");
}
