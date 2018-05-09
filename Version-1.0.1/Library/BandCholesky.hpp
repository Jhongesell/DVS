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


#ifndef __BandCholesky__
#define __BandCholesky__

#include "Definiciones.hpp"
#include "Solvable.hpp"
#include "MatrizDispersa.hpp"
#include "ErrorControl.hpp"



/* Here a symmetric banded matrix aij is represented as a set of bw + 1 columns by n
   rows. The diagonal aii is the column A[bw][i], and aij is A[bw + j - i][i] for j <= i
   where i - j <= bw. The matrix AK[bw][n] is overwritten by the factorLU() method.
*/

class BandCholesky : public Solvable
{

protected:

   int n;
   int bw;
   ldouble **AK;
   /// Control de errores
   ErrorControl ce;


   void factorLU(void);


public:


   /* Convert square matrix A into a Cholseky banded matrix and then
      do a standard BandCholesky. Leaves bandwidth in global variable bw  */
   BandCholesky(int n, MatrizDispersa *A)
   {
      AK = NULL;
      name = "BandCholesky";
      convertBand(n, A);
      this->n = n;
      factorLU();
   }

   ~BandCholesky()
   {
      clean();
   }

   void clean(void)
   {
      if (AK == NULL) return;
      int i;
      for (i = 0; i < (bw + 1); i++)
      {
         delete []AK[i];
         AK[i] = NULL;
      }
      delete []AK;
      AK = NULL;
   }


   /* returns a banded Cholesky matrix from a square symmetrix matrix */
   void convertBand(int n, ldouble **A);

   /* returns a banded Cholesky matrix from a square symmetrix matrix */
   void convertBand(int n, MatrizDispersa *A);


   void solve(ldouble *x, ldouble *y);

   void print(void);

   inline int getIter(void)
   {
      return 0;
   }

};

#endif
