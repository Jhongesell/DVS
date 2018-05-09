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



#include <stdio.h>
#include "Definiciones.hpp"
#include "PropDef.hpp"
#include "LM1MPI.hpp"
#include "PLM1MPI.hpp"
#include "LM2MPI.hpp"
#include "PLM2MPI.hpp"
#include "MF1MPI.hpp"
#include "PMF1MPI.hpp"
#include "MF2MPI.hpp"
#include "PMF2MPI.hpp"
#include "ErrorControl.hpp"



// Programa Maestro-Esclavo
int main(int argc, char *argv[])
{
   /// Control de errores
   ErrorControl ce("", "main");

   // Variables para MPI
   int MPI_id, MPI_np;
   MPI::Init(argc, argv);
   MPI_id = MPI::COMM_WORLD.Get_rank();
   MPI_np = MPI::COMM_WORLD.Get_size();


   // Revisa que pueda arrancar el esquema M-E
   if (MPI_np < 2)
   {
      printf("Se necesitan almenos dos procesadores para el esquema M-E\n");
      return 1;
   }

   // Control de parametros pasados de la linea de comandos
   PropDef *props = new (nothrow) PropDef();
   if (props == NULL) ce.memoryError("props");
   props->parse(argc, argv);


   ldouble a[3];
   ldouble b[3];
   ldouble c;
   char *method;

   method = props->getString("m", "CMG");

   int nDim = props->getInt("d", 0);
   if (nDim < 2 || nDim > 3) ce.fatalError(1, "No se especifico la dimension del problema");

   a[0] = props->getDouble("ax", 0.0);
   a[1] = props->getDouble("ay", 0.0);
   a[2] = props->getDouble("az", 0.0);
   b[0] = props->getDouble("bx", 0.0);
   b[1] = props->getDouble("by", 0.0);
   b[2] = props->getDouble("bz", 0.0);
   c = props->getDouble("c", 0.0);
   EllipOp *op = new (nothrow) EllipOp(nDim, a, b, c);
   if (op == NULL) ce.memoryError("op");

   // Ejecuta el metodo correspondiente
   if (strcmp(method, "LM1") == 0)
   {
      LM1MPI *dualp = new (nothrow) LM1MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "LM2") == 0)
   {
      LM2MPI *dualp = new (nothrow) LM2MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "PLM1") == 0)
   {
      PLM1MPI *dualp = new (nothrow) PLM1MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "PLM2") == 0)
   {
      PLM2MPI *dualp = new (nothrow) PLM2MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "MF1") == 0)
   {
      MF1MPI *dualp = new (nothrow) MF1MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "MF2") == 0)
   {
      MF2MPI *dualp = new (nothrow) MF2MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "PMF1") == 0)
   {
      PMF1MPI *dualp = new (nothrow) PMF1MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else if (strcmp(method, "PMF2") == 0)
   {
      PMF2MPI *dualp = new (nothrow) PMF2MPI(MPI_id, MPI_np, *props, *op);
      if (dualp == NULL) ce.memoryError("dualp");
      dualp->solvePar();
      dualp->clean();
      delete dualp;
   }
   else  printf("\nNo valid method specified");

   MPI::Finalize();
   printf("\n\nEnd Program\n\n");
   fflush(stdout);

   delete []method;
   delete op;
   delete props;


   return 0;
}

