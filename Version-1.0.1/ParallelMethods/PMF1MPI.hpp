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



#ifndef __PMF1MPI__
#define __PMF1MPI__


#include "DPMainMPI.hpp"
#include "InterchangeMPI.hpp"
#include "PMF1.hpp"


/// Clase para definir el metodo PMF-1 de DVS-DDM
/** Clase para definir el metodo PMF-1 de DVS-DDM en paralelo
*/
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.0
    @bug No hay errores conocidos
*/
class PMF1MPI: public DPMainMPI, public  PMF1
{

public:

   /// Constructor de la clase
   PMF1MPI(int id, int np, PropDef &props, EllipOp &op) : DPMainMPI(id, np, props, op), PMF1(props, op)
   {
   }

   /// Inicializa InterchangeMPI en lugar de Interchange
   void iniInterchage(void)
   {
      DPMainMPI::ce.nameClassFunct("PMF1MPI", "iniInterchage");

      EsquemaMEMPI *me = this;
      inter = new InterchangeMPI(*props, *me);
      if (inter == NULL) DPMainMPI::ce.memoryError("inter");
   }

   void clean(void)
   {
      //~ solver->clean();
      delete solver;
      solver = NULL;
   }

   /// Sobrecarga del la aplicacion
   void solvePar(void)
   {
      if (id == 0)
      {
         // Manda el indice al cual hay que crear la geometria
         for (int i = 0; i < DPMainMPI::nOmega; i++)
         {
            reparteCargaTrabajo(xnp, indl, i);
            mss[0] = indl;
            mss[1] = i;
            //~ printf("\nEnvio %d %d %d \n", i, indl, xnp);
            //~ fflush(stdout);
            MPI::COMM_WORLD.Send(&mss, 2, MPI::INT, xnp, 1);
         }

         PMF1MPI::solve();
         fflush(stdout);

         // Finalizacion de las tareas
         mss[0] = 0;
         mss[1] = 0;
         mss[2] = 0;
         for (int i = 1; i < np; i++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, i, 1);
      }
      else
      {
         if (nta) Esclavo();
         else
         {
            int sw = 1;
            while (sw)
            {
               MPI::COMM_WORLD.Recv(&msa, 10, MPI::INT, 0, 1);
               sw = msa[2];
            }
         }
      }
   }



};


#endif


