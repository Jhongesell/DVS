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


#ifndef __ICGM__
#define __ICGM__

#include "Definiciones.hpp"
#include "CGM.hpp"
#include "MatrizDispersa.hpp"



/// Clase para implementar CGM con matrices bandadas o dispersas
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ICGM : public CGM, public  MultOp, public DotProd
{

private:

   /// Multiplica Au=v
   MatrizDispersa *M;

   /// Variables temporales
   ldouble val;
   int i;

   /// Producto punto
   inline ldouble dot(ldouble *u, ldouble *v)
   {
      val = 0.0;
      for (i = 0; i < n; i++) val += u[i] * v[i];
      return val;
   }

   /// Multiplica Au=v
   inline void multOp(ldouble *u, ldouble *v)
   {
      M->multiplica(u, v);
   }


public:

   /// Contructor de la clase
   ICGM(int n, MatrizDispersa *M, ldouble eps, int iter)
   {
      CGM::n = n;
      CGM::eps = eps;
      CGM::nMaxIter = iter;
      this->M = M;
      A = (MultOp*) this;
      dotP = (DotProd*) this;
      inicializa();
   }

   /// Destructor de la clase
   ~ICGM()
   {
      clean();
   }

   void clean(void)
   {
      delete M;
      M = NULL;
   }

   /// vector size
   inline int getSize(void)
   {
      return n;
   }

};

#endif

