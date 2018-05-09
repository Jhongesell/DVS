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



#ifndef  __DualPrimal__
#define __DualPrimal__

#include <vector>
#include "Definiciones.hpp"
#include "BdNode.hpp"
#include "Solvable.hpp"
#include "Interchange.hpp"
#include "ErrorControl.hpp"



class DualPrimal
{
protected:

   int nPrimal;
   int nDual;
   Interchange *inter;
   vector<Solvable*> SP;                   // SP[0] = internal node inverse, SP[1] = complete inverse

   ldouble *XP;                             // Solution of SP[i].solve(XP, YP)
   ldouble *YP;                             // RHS of SP[i].solve(XP, YP)

   int nD;                                 // size of bdDuals
   int nP;                                 // size of bdPrimals
   /// Control de errores
   ErrorControl ce;



public:

   DualPrimal(Interchange &inter);


   ~DualPrimal(void)
   {
      size_t i;
      for (i = 0; i < SP.size(); i++)
      {
         Solvable *abc = SP[i];
         if (abc != NULL) abc->clean();
         delete abc;
         abc = NULL;
      }

      delete []XP;
      delete []YP;
   }


   void a(ldouble *u, ldouble *v) ;

   void calcValues(ldouble *u);

   void fromSubdomains(int sc, ldouble *u);

   void genMats(void);

   inline int getNDual(void)
   {
      return nDual;
   }


   void j(ldouble *u, ldouble *v);

   void multS(ldouble *u, ldouble *v);


   void solveAPP(int sp, int sc1, int sc2, int sc3);

   void solveS(ldouble *u, ldouble *v);

   void toSubdomains(int sc, ldouble *u);


};

#endif



