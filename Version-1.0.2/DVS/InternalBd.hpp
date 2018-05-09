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


#ifndef __InternalBd__
#define __InternalBd__



/* This class represents an internal boundary element. A list of these
   elements/subdomain is generated by the Subdomain class on a call to
   the method getInternalBd()
*/
#include <math.h>
#include <string.h>
#include "Definiciones.hpp"
#include "ErrorControl.hpp"


class InternalBd
{

private:

   int subd;          // subdomain number
   int node;          // node number (of the subdomain)
   int bd;            // sequential boundary node number
   int dp;            // 0 if the node is dual, -1 if the node is primal
   int nDim;          // problem space dimension
   ldouble *coord;     // [0:nDim - 1] array of node spatial coordinates
   //~ /// Control de errores
   //~ ErrorControl ce;



public:

   InternalBd(void)
   {
      coord = NULL;
   }

   InternalBd(int s, int n, int b, int i, int d, ldouble *cor)
   {
      subd = s;
      node = n;
      bd = b;
      dp = i;
      nDim = d;
      coord = cor;
   }

   ~InternalBd(void)
   {
      delete []coord;
      coord = NULL;
   }


   bool equals(InternalBd *x)
   {
      int i;
      for (i = 0; i < nDim; i++)
      {
         //~ printf("\n%d ==> %f %f",i, coord[i], x->rcoord(i));
         if (fabs(coord[i] - x->rcoord(i)) > EPS_EQUAL) return false;
      }
      return true;
   }

   int compareTo(InternalBd *a)
   {
//~ printf("\n %d %d %d %d",subd, node, dp,nDim);
      int r = 0, i;
      int n = nDim - 1;
      for (i = n; i >= 0; i--)
      {
//~ printf("\n %f %f",coord[i],a->rcoord(i));
         if (coord[i] < a->rcoord(i) - EPS_EQUAL)
         {
            r = -1;
            break;
         }
         else if (a->rcoord(i) < coord[i]  - EPS_EQUAL)
         {
            r = 1;
            break;
         }
      }
//~ printf("\n%d",r);
      if (r != 0) return r;
      if (subd < a->rsubd()) return -1;
      if (subd > a->rsubd()) return 1;
      return 0;
   }


   /*
      char* toString(void)
      {
         ce.nameClassFunct("InternalBd","toString");

         char *r = new (nothrow) char[300];
         if (r == NULL) ce.memoryError("r");
         char *p = new (nothrow) char[50];
         if (p == NULL) ce.memoryError("p");

         sprintf(r,"[ %d %d %d %d %d (",subd, node, bd, dp, nDim);
         //~ String s = "[" + subd + " " + node + " " + bd + " " + dp + " " + nDim + " (";
         for (int i = 0; i < nDim; i++)
         {
   #ifdef __Double__
            sprintf(p,"%f",coord[i]);
   #else
            sprintf(p,"%Lf",coord[i]);
   #endif
            strcat(r,p);
            if (i == nDim - 1) strcat(r,")");
            else strcat(r," ");
         }
         strcat(r,"]");
         delete []p;
         return r;
      }
   */

   // Se crean para poder acceder a datos privados de la clase
   inline ldouble rcoord(int i)
   {
      return coord[i];
   }

   inline int rsubd(void)
   {
      return subd;
   }

   inline int rnode(void)
   {
      return node;
   }

   inline int rdp(void)
   {
      return dp;
   }

   inline int rbd(void)
   {
      return bd;
   }


   void getval(int &s, int &n, int &b, int &i, int &d, ldouble *c)
   {
      s = subd;
      n = node;
      b = bd;
      i = dp;
      d = nDim;
      int j;
      for (j = 0; j < nDim; j++) c[j] = coord[j];
   }

};

#endif