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



#ifndef __BandSolve__
#define __BandSolve__

#include "Definiciones.hpp"
#include "Solvable.hpp"
#include "MatrizDispersa.hpp"
#include "ErrorControl.hpp"



/* Here a banded matrix aj is represented as a set of (2*bw + 1) columns by n
   rows. The diagonal aii is the column A[bw][i], and aij is A[bw + j - i][i]
   where |i - j| <= bw.
*/

class BandSolve : public Solvable
{

private:

   int bw;          // bandwidth
   int n;           // matrix length
   ldouble **AK;      // AK[2*bw + 1][n]


protected:

   // Control de errores
   ErrorControl ce;


   void factorLU(void);


public:

   BandSolve(void)
   {
      name = "BandSolve";
      bw = n = 0;
      AK = NULL;
   }

   BandSolve(int n, ldouble **A);

   BandSolve(int n, MatrizDispersa *A);

   ~BandSolve()
   {
      clean();
   }

   void clean(void)
   {
      if (AK == NULL) return;
      int i;
      for (i = 0; i < (2 * bw + 1); i++)
      {
         delete []AK[i];
         AK[i] = NULL;
      }
      delete []AK;
      AK = NULL;
   }

   void solve(ldouble *x, ldouble *y);

   void convertBand(int n, ldouble **A);

   void convertBand(int n, MatrizDispersa *A);

   void print(void);

   inline int getIter(void)
   {
      return 0;
   }
};

#endif
