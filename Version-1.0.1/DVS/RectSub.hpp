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



#ifndef __RectSub__
#define __RectSub__

#include <vector>
#include "Definiciones.hpp"
#include "InternalBd.hpp"
#include "EllipOp.hpp"
#include "Solvable.hpp"
#include "Primal.hpp"
#include "MatrizDispersa.hpp"
#include "ErrorControl.hpp"

using namespace std;


class RectSub //: public Subdomain
{
protected:

   static const int KNOWN = 1, INTERIOR = 2, INTBD = 4, VERTEX = 8, EDGE = 16, FACE = 32, PRIMAL = 64, DUAL = 128;


   ldouble Ca[3];
   ldouble Cb[3];
   ldouble  Cc;

   int id;                     // Subdomain number
   int nDim;                   // space dimension (1, 2 or 3)
   ldouble **domain;          // [nDim][2]  min, max in each dimension
   EllipOp *op;
   int *mesh;
   int *N;
   int *M;                    // mesh size [number of subdivisions per axis]
   int *M1;
   int *coord;                // integer coordinates of node number
   int *coordN;               // integer coordinates of node number (coarse grid)
   ldouble *h;                 // mesh size for each dimension
   ldouble hfac;                // element measure h[0]*h[1]*...
   ldouble **scr;             // scratch arrays for (external) calculations
   int np;                     // total node size (Mx +1) or (Mx + 1)(My + 1) or (Mx + 1)*(My + 1)*(Mz + 1)
   int *ntype;                // Describes type of each node
   vector<Solvable*> inv;             // Inversea of interiors and complete (minus primals) matrices
   ldouble *coef;              // coefficient vector (for the operator)
   int *bdMap;                // Mapping between serial bd# and local node#
   int *mapInt;               // equation mapping for interior matrix
   int *mapFull;              // equation mapping for full (minus primals) matrix
   int nInt;                   // number of interior nodes (for mapInt)
   int nFull;                  // number of interior + dual nodes
   ldouble *X;                 // Used for computing inverses (result)
   ldouble *Y;                 // used for computing inverses r.h.s
   bool bFloat;             // true if subdomain has no known values (ie no external boundary)
   bool bsym;               // true if operator is symmetric
   ldouble *x;                 // (x,y,z) coordinates for a node number (computed by getCoord())
   int nBd;
   /// Control de errores
   ErrorControl ce;

public:

   RectSub(int id, int nDim, int *mesh, ldouble **dom, EllipOp &op, Primal &primal);


   ~RectSub(void)
   {
      int i;

      for (i = 0; i < 4; i++)
      {
         delete []scr[i];
         scr[i] = NULL;
      }
      delete []scr;
      scr = NULL;

      for (i = 0; i < nDim; i++)
      {
         delete []domain[i];
         domain[i] = NULL;
      }
      delete []domain;
      domain = NULL;

      delete []mapInt;
      delete []mapFull;
      delete []X;
      delete []Y;
      delete []x;
      delete []N;
      delete []M;
      delete []M1;
      delete []coord;
      delete []coordN;
      delete []h;
      delete []ntype;
      delete []coef;
      delete []bdMap;

      size_t j;
      for (j = 0; j < inv.size(); j++)
      {
         Solvable *abc = inv[j];
         abc->clean();
         delete abc;
         abc = NULL;
      }
   }

   int addProjNs(ldouble **A, int *map, ldouble fac);

   int addProjNs(MatrizDispersa *A, int *map, ldouble fac);

   void clear(int s);

   void diff(int sc3, int sc1, int sc2);

   inline void diffValues(int sc, ldouble *u)
   {
      for (int i = 0; i < nBd; i++) u[i] -= scr[sc][bdMap[i]];
   }


   void genCoef(EllipOp &op);

   void genCoefVar(int ren);

   void genInv(int type);

   Solvable *genInverse(int *map, ldouble fac);

   void genNcoord(int n, int *coord, int *N);

   void genNtype(Primal &primal);

   inline int getBdSize()
   {
      return nBd;
   }

   void getCoord(int m, ldouble *x);

   inline void getCoordNode(int n, ldouble *x)
   {
      getCoord(bdMap[n], x);
   }

   vector<InternalBd*> getInternalBd(void);

   inline vector<Solvable*> getInv(void)
   {
      return inv;
   }

   inline int *getNtype(void)
   {
      return ntype;
   }

   // Se agrega para poder actualizar Ntype en paralelo (ya que no se pasan datos por referencia)
   inline void setNtype(int *arr)
   {
      for (int i = 0; i < np; i++) ntype[i] = arr[i];
   }

   void getPrimals(int sc, ldouble *u);

   inline ldouble getValue(int sc, int n)
   {
      return scr[sc][bdMap[n]];
   }

   void getValues(int sc, ldouble *u);

   void inverse(int sp, int sc1, int sc2);

   bool isKnown(int *coord);

   bool isInterior(int *coord);

   bool isIntBd(int *coord);

   int nodeType(int *coord);

   inline bool isDual(int i)
   {
      return (ntype[i] & DUAL) != 0;
   }

   inline bool isFloat(void)
   {
      return bFloat;
   }


   inline bool isInterior(int i)
   {
      return (ntype[i] & INTERIOR) != 0;
   }

   inline bool isKnown(int i)
   {
      return (ntype[i] & KNOWN) != 0;
   }

   inline bool isPrimal(int i)
   {
      return (ntype[i] & PRIMAL) != 0;
   }

   inline bool isVertex(int i)
   {
      return (ntype[i] & VERTEX) != 0;
   }

   void knownValues(int s1);

   void multOp(int s1, int s2);

   void printMat(const char *s, ldouble **A, int tm);

   void printMult(void);

   void rhs(int sc);

   void setPrimals(int sc, ldouble *u);

   inline void setValue(int sc, int n, ldouble val)
   {
      scr[sc][bdMap[n]] = val;
   }

   void setValues(int sc, ldouble *u);

   void print(const char *s, int sc);

   void print(int sc);

   inline int getNP(void)
   {
      return np;
   }
};


#endif

