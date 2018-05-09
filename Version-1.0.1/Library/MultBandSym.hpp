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



#ifndef __MultBandSym__
#define __MultBandSym__


#include "Definiciones.hpp"
#include "MultOp.hpp"


/* Multiplication algorithm for banded symmetric matrices. aij = AK[bw + j - i][i]  j<= i
*/
class MultBandSym : public MultOp
{
public:

   int n;              // Matrix size
   int bw;             //
   ldouble **AK;        // AK[bw + 1][n];   lower band

   MultBandSym(int n, int bw, ldouble ** AK)
   {
      this->AK = AK;
      this->n = n;
      this->bw = bw;
   }

   void multOp(ldouble *x, ldouble *y);

   inline int getSize(void)
   {
      return n;
   }

};

#endif
