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



#ifndef __IDQGMRES__
#define __IDQGMRES__


#include "Definiciones.hpp"
#include "DQGMRES.hpp"
#include "MatrizDispersa.hpp"



/// Clase para implementar DQGMRES con matrices bandadas o dispersas
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.1
    @bug No hay errores conocidos
*/
class IDQGMRES : public DQGMRES, public MultOp
{

private:

   /// Matriz Bandada o Dispersa
   MatrizDispersa *M;

   /// Multiplica Au=v
   inline void multOp(ldouble *u, ldouble *v)
   {
      M->multiplica(u, v);
   }


public:

   /// Constructor de la clase
   IDQGMRES(int n, MatrizDispersa *M, int k, double eps, int iter)
   {

      this->M = M;
      mult = (MultOp*) this;
      DQGMRES::n = n;
      DQGMRES::k = k;
      DQGMRES::eps = eps;
      DQGMRES::nMaxIter = iter;
      inicializa();
   }

   // Destructor de la clase
   ~IDQGMRES()
   {
      clean();
   }

   void clean(void)
   {
      delete M;
      M = NULL;
   }

   /// vector size
   int getSize(void)
   {
      return n;
   }

};

#endif

