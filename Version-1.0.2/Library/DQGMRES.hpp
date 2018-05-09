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



#ifndef __DQGMRES__
#define __DQGMRES__

#include <math.h>
#include "Definiciones.hpp"
#include "Solvable.hpp"
#include "MultOp.hpp"
#include "ErrorControl.hpp"





class DQGMRES : public Solvable
{
protected:

   int n;             // vector size
   int k;             // # vectors to save
   int k1;
   int maxIter;       // maximum number of interations
   int nIter;         // current iteration
   MultOp *mult;       // Operator to solve A*x = b
   ldouble gm;         // gamma(m)
   ldouble gm1;        // gamma(m + 1)
   ldouble **p;      // [k + 1][n]
   ldouble **cs;     // [k + 1][2] cosine, sine array for rotation
   ldouble **h;      // Hessenberg
   ldouble **q;      // quasi orthogonal vectors [k + 1][n]
   ldouble *v;        // scratch
   ldouble eps;
   int nMaxIter;
   // Control de errores
   ErrorControl ce;


public:

   /* n = vector size, k = number of previous vectors to be saved & orthogonalized  */
   DQGMRES(MultOp &mult, int k, ldouble eps)
   {
      name = "DQGMRES";
      this->eps = eps;
      this->nMaxIter = NMAXITER;
      this->mult = &mult;
      this->n = mult.getSize();
      this->k = k;
      inicializa();
   }

   DQGMRES(void)
   {
      this->eps = EPSILON;
      this->nMaxIter = NMAXITER;
   }

   ~DQGMRES(void)
   {
      clean();
   }

   void clean(void)
   {
      int i;
      if (p)
      {
         for (i = 0; i < k1; i++)
         {
            delete []p[i];
            p[i] = NULL;
         }
         delete []p;
         p = NULL;
      }

      if (q)
      {
         for (i = 0; i < (k1 + 1); i++)
         {
            delete []q[i];
            q[i] = NULL;
         }
         delete []q;
         q = NULL;
      }

      if (cs)
      {
         for (i = 0; i < k1; i++)
         {
            delete []cs[i];
            cs[i] = NULL;
         }
         delete []cs;
         cs = NULL;
      }

      if (h)
      {
         for (i = 0; i < (k1 + 1); i++)
         {
            delete []h[i];
            h[i] = NULL;
         }
         delete []h;
         h = NULL;
      }

      delete []v;
      v = NULL;
   }


   void inicializa(void);

   void applyOmega(int m);

   void solve(ldouble *x, ldouble *b);

   inline int getIter(void)
   {
      return nIter;
   }

   inline void setMaxIter(int nmi)
   {
      nMaxIter = nmi;
   }

   inline void setEpsilon(ldouble ep)
   {
      eps = ep;
   }

};

#endif
