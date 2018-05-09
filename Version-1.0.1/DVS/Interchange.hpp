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



#ifndef __Interchange__
#define __Interchange__

#include "Definiciones.hpp"
#include "BdNode.hpp"
#include "Solvable.hpp"
#include "PropDef.hpp"
#include "RectSub.hpp"
#include "FunctionV1.hpp"
#include "Constant.hpp"
#include "EllipOp.hpp"
#include "Primal.hpp"
#include "LookUpFunction.hpp"
#include "VertPrimal.hpp"
#include "VertEdgePrimal.hpp"
#include "AllPrimal.hpp"
#include "NoPrimal.hpp"
#include "CreateBdNodes.hpp"
#include "ErrorControl.hpp"




/* Controls all communication and data exchange between master and subdomains */
class Interchange
{
protected:

   vector<RectSub*> omegas;
   ldouble** bdValues;

   int nOmega;
   int nDim;
   PropDef *props;
   EllipOp *op;

   FunctionV1 *zero, *one;
   FunctionV1 *f, *g;
   char *sf, *sg;
   ldouble fc, gc;
   int *mesh;
   char *prim;
   char *method;
   int swprint;
   ldouble Ax, Ay, Az;

   ldouble **domain;
   ldouble ax, ay, az;
   ldouble c;
   ldouble bx, by, bz;
   Primal *primal;
   /// Control de errores
   ErrorControl ce;

   // generate the geometry
   void  genGeom(void);


public:

   CreateBdNodes *bds;


   /// Constructor
   Interchange(PropDef &props);

   /// Constructor
   Interchange(void)
   {
      bdValues = NULL;
      props = NULL;
      op = NULL;
      zero = NULL;
      one = NULL;
      f = NULL;
      g = NULL;
      sf = NULL;
      sg = NULL;

      mesh = NULL;
      prim = NULL;
      method = NULL;

      domain = NULL;
      primal = NULL;
      bds = NULL;
   }

   /// Destructor
   virtual ~Interchange()
   {
      int i;

      delete bds;


      for (i = 0; i < nOmega; i++)
      {
         delete []bdValues[i];
      }
      delete []bdValues;
      bdValues  = NULL;

      delete []mesh;
      mesh = NULL;

      for (i = 0; i < 3; i++) delete [] domain[i];
      delete []domain;

      delete zero;
      zero = NULL;
      delete one;
      one = NULL;

      delete f;
      f = NULL;
      delete g;
      g = NULL;

      delete []sf;
      delete []sg;
      delete []prim;
      delete []method;
      delete op;
      delete primal;

      size_t j;
      for (j = 0; j < omegas.size(); j++)
      {
         RectSub *abc = omegas[j];
         delete abc;
         abc = NULL;
      }
      omegas.clear();
   }


   void initialize(int nOmega);


   /* Computes APP-1 where sp = 0 for APP with interior nodes
                           sp = 1 for APP with interior + dual nodes
      This method is parallelizable to a point.
   */
   Solvable *calcSP(int sp);

   inline int getND(void)
   {
      return bds->nD;
   }

   inline int getNP(void)
   {
      return bds->nP;
   }

   inline int getNPrimals(void)
   {
      return bds->nPrimal;
   }

   inline int getNDuals(void)
   {
      return bds->nDual;
   }

   inline int getnOmega(void)
   {
      return nOmega;
   }

   inline int getnDim(void)
   {
      return nDim;
   }


   inline ldouble rbdValues(int i, int j)
   {
      return bdValues[i][j];
   }

   inline void sbdValues(int i, int j, ldouble v)
   {
      bdValues[i][j] = v;
   }

   void pbdValues(void)
   {
      int i, j;
      for (i = 0; i < nOmega; i++)
      {
         printf("\n");
#ifdef __Double__
         for (j = 0; j < bds->maxBd; j++) printf(" %e ", bdValues[i][j]);
#else
         for (j = 0; j < bds->maxBd; j++) printf(" %Le ", bdValues[i][j]);
#endif
      }
   }

   /// Clear scr[sc][] en e subdomains
   void clear(int e, int sc)
   {
      omegas[e]->clear(sc);
   }

   void setValue(int e, int sc, int n, ldouble val)
   {
      omegas[e]->setValue(sc, n, val);
   }

   /// scr[sc2][] = A(sp)-1(scr[sc1][])
   void inverse(int e, int sp, int sc1, int sc2)
   {
      omegas[e]->inverse(sp, sc1, sc2);
   }

   /// scr[s2][] = A(scr[sc1][])
   void multOp(int e, int sc1, int sc2)
   {
      omegas[e]->multOp(sc1, sc2);
   }

   virtual void calcula(int e, int node, int sp)
   {
      clear(e, 0);
      setValue(e, 0, node, 1.0);
      multOp(e, 0, 1);
      inverse(e, sp, 1, 2);
      multOp(e, 2, 0);
   }

   /// Clear scr[sc][] in all subdomains
   virtual void clear(int sc)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->clear(sc);
   }

   virtual ldouble getValue(int e, int scr, int node)
   {
      return omegas[e]->getValue(scr, node);
   }

   virtual ldouble getValue(int e, int scr1, int scr2, int node)
   {
      return (omegas[e]->getValue(scr1, node) - omegas[e]->getValue(scr2, node));
   }

   /// scr[sc3][] = scr[sc1][] - scr[sc2][] in all subdomains
   virtual void diff(int sc3, int sc1, int sc2)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->diff(sc3, sc1, sc2);
   }

   /// scr[sc2][] = A(sp)-1(scr[sc1][])
   virtual void inverse(int sp, int sc1, int sc2)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->inverse(sp, sc1, sc2);
   }

   /// scr[sc][] = Dirichlet boundary values of all subdomains
   virtual void knownValues(int sc)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->knownValues(sc);
   }

   /// scr[s2][] = A(scr[sc1][])
   virtual void multOp(int sc1, int sc2)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->multOp(sc1, sc2);
   }

   /// scr[sc][] = initial right-hand-side (all subdomains)
   virtual void rhs(int sc)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->rhs(sc);
   }


   virtual void genInv(int e, int type)
   {
      omegas[e]->genInv(type);
   }

   virtual void getCoordNode(int e, int n, ldouble *x)
   {
      omegas[e]->getCoordNode(n, x);
   }

   virtual void print(const char *s, int sc)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->print(s, sc);
   }

   virtual void print(int sc)
   {
      for (int i = 0; i < nOmega; i++) omegas[i]->print(sc);
   }

   virtual int getMaxBdSize(void)
   {
      int n, e, maxBd = 0;
      for (e = 0; e < nOmega; e++)
      {
         n = omegas[e]->getBdSize();
         if (maxBd < n) maxBd = n;
      }
      return maxBd;
   }

   virtual int *getNtype(int e)
   {
      return omegas[e]->getNtype();
   }

   // Se agrega para poder actualizar Ntype en paralelo (ya que no se pasan datos por referencia)
   virtual void setNtype(int e, int *arr)
   { }



   /// bdValues[][] -= scr[sc][] in all subdomains
   virtual void diffValues(int sc)
   {
      for (int s = 0; s < nOmega; s++) omegas[s]->diffValues(sc, bdValues[s]);
   }

   /// bdValues[][] = scr[sc][] from all subdomains
   virtual void fromSubdomains(int sc)
   {
      for (int s = 0; s < nOmega; s++) omegas[s]->getValues(sc, bdValues[s]);
   }

   /// bdValues[][] (primals only) = scr[sc][] (primals)
   virtual void getPrimals(int sc)
   {
      for (int s = 0; s < nOmega; s++) omegas[s]->getPrimals(sc, bdValues[s]);
   }

   /// scr[sc][] = bdValues  all subdomains
   virtual void setPrimals(int sc)
   {
      for (int s = 0; s < nOmega; s++) omegas[s]->setPrimals(sc, bdValues[s]);
   }

   /// scr[sc][] = bdValues[][]  all subdomains
   virtual void toSubdomains(int sc)
   {
      for (int s = 0; s < nOmega; s++) omegas[s]->setValues(sc, bdValues[s]);
   }


};

#endif
