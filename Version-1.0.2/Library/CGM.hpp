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


#ifndef __CGM__
#define __CGM__

#include "Definiciones.hpp"
#include "Solvable.hpp"
#include "MultOp.hpp"
#include "DotProd.hpp"
#include "ErrorControl.hpp"



class CGM : public Solvable
{

protected:
   // size of vector
   int n;
   // scratch vectors
   ldouble *r, *p, *v;
   // scratch scalars
   ldouble alpha, beta, gamma, lambda;
   // allowed error
   ldouble eps;

   ldouble mu;
   // Operator
   MultOp *A;
   // Inner Product
   DotProd *dotP;
   // number of iterations
   int nIter;
   // Numero maximo de iteraciones
   int nMaxIter;
   // Control de errores
   ErrorControl ce;

   ldouble norm(ldouble *x);

public:

   CGM(void)
   {
      nMaxIter = NMAXITER;
      this->eps = EPSILON;
      r = NULL;
      p = NULL;
      v = NULL;
   }


   CGM(MultOp &A, DotProd &dotP, ldouble eps)
   {
      ce.nameClassFunct("CGM", "CGM");
      r = NULL;
      p = NULL;
      v = NULL;

      name = "CGM";
      nMaxIter = NMAXITER;
      this->A = &A;
      this->dotP = &dotP;
      this->eps = eps;
      n = A.getSize() * DIM_VECTOR;
      inicializa();
   }


   ~CGM()
   {
      clean();
   }

   void clean(void)
   {
      delete []r;
      r = NULL;
      delete []p;
      p = NULL;
      delete []v;
      v = NULL;
   }


   void inicializa(void)
   {
      r = new (nothrow) ldouble[n];
      if (r == NULL) ce.memoryError("r");
      p = new (nothrow) ldouble[n];
      if (p == NULL) ce.memoryError("p");
      v = new (nothrow) ldouble[n];
      if (v == NULL) ce.memoryError("v");
      for (int i = 0; i < n; i++) r[i] = p[i] = v[i] = 0.0;
   }


   void solve(ldouble *u, ldouble *b);


   int getIter(void)
   {
      return nIter;
   }


   void setMaxIter(int nmi)
   {
      nMaxIter = nmi;
   }

   void setEpsilon(ldouble ep)
   {
      eps = ep;
   }

};

#endif
