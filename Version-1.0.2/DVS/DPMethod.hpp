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


#ifndef __DPMethod__
#define __DPMethod__


#include "Definiciones.hpp"
#include "MultOp.hpp"
#include "DotProd.hpp"
#include "PropDef.hpp"
#include "DualPrimal.hpp"
#include "Interchange.hpp"
#include "ErrorControl.hpp"



class DPMethod : public  MultOp, public DotProd
{

protected:

   PropDef *props;                        // Property list passed from Problem
   int printv;                            // Print (0 = no print, 1 = dual values, 2 = all values, 3 all values)
   ldouble epsilon;                        // Allowed error for residual
   int nDual;                             // Number of dual nodes
   int nOmega;                            // Number of Subdomains
   int nDim;                              // Space dimension
   DualPrimal *dualp;                     // Dual primal calculations
   Interchange *inter;

   ldouble *u;                             // Equation solution on the interior boundary
   ldouble *rhss;                          // Construct r.h.s of equation
   ldouble *scr;                           // Scratch array
   Solvable *solver;                      // Object that solves the equation

   time_t time0;                          // Start of processessing
   time_t time1;                          // End of local inverse calculation
   time_t time2;                          // End of primal matrix calculation
   time_t time3;                          // End of iterations

   /// Control de errores
   ErrorControl ce;


   /// Inicializa los subdominios
   virtual void iniInterchage(void)
   {
      ce.nameClassFunct("DPMethod", "iniInterchage");

      inter = new (nothrow) Interchange(*props);
      if (inter == NULL) ce.memoryError("inter");
   }


public:

   // Constructor
   DPMethod(PropDef &props)
   {
      scr = NULL;
      inter = NULL;
      epsilon = EPSILON;
      u = NULL;
      rhss = NULL;
      dualp = NULL;
      solver = NULL;
      this->props = &props;
   }

   virtual ~DPMethod()
   {
      delete []u;
      u = NULL;

      delete []rhss;
      rhss = NULL;

      delete dualp;
      dualp = NULL;
   }

   void initialize(void);

   virtual void clean(void) = 0;

   void genInverse(int type);

   inline int getSize(void)
   {
      return nDual;
   }

   void print(ldouble *u);

   void printTime(void);

   const char *prCoord(ldouble *x);

   virtual void rhs(void) = 0;

   virtual void solve(void) = 0;

   double analyticSolution(double *x);

   void conditionalNumber(bool symetric);
};

#endif

