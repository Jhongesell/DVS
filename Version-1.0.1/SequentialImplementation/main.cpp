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


#include <string.h>
#include <time.h>
#include "Definiciones.hpp"
#include "PropDef.hpp"
#include "MF1.hpp"
#include "PMF1.hpp"
#include "MF2.hpp"
#include "PMF2.hpp"
#include "LM1.hpp"
#include "PLM1.hpp"
#include "LM2.hpp"
#include "PLM2.hpp"
#include "ErrorControl.hpp"


int main(int argc, char *args[])
{
   /// Control de errores
   ErrorControl ce("", "main");

   time_t t1, t2;
   time(&t1);

   PropDef *props = new (nothrow) PropDef();
   if (props == NULL) ce.memoryError("props");
   props->parse(argc, args);

   ldouble a[3];
   ldouble b[3];
   ldouble c;
   char *method;

   DPMethod *dualp;
   method = props->getString("m", "MF1");

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

   dualp = NULL;
   if (strcmp(method, "MF1") == 0)         dualp = new (nothrow) MF1(*props, *op);
   else if (strcmp(method, "PMF1") == 0)   dualp = new (nothrow) PMF1(*props, *op);
   else if (strcmp(method, "MF2") == 0)    dualp = new (nothrow) MF2(*props, *op);
   else if (strcmp(method, "PMF2") == 0)   dualp = new (nothrow) PMF2(*props, *op);
   else if (strcmp(method, "LM1") == 0)    dualp = new (nothrow) LM1(*props, *op);
   else if (strcmp(method, "PLM1") == 0)   dualp = new (nothrow) PLM1(*props, *op);
   else if (strcmp(method, "LM2") == 0)    dualp = new (nothrow) LM2(*props, *op);
   else if (strcmp(method, "PLM2") == 0)   dualp = new (nothrow) PLM2(*props, *op);
   else ce.fatalError(1, "No valid method specified");
   if (dualp) dualp->solve();


   time(&t2);
   printf("\nTiempo Calculo: %f\n\n\n", difftime(t2, t1));

   dualp->clean();
   delete dualp;
   delete []method;
   delete props;
   delete op;


   return 0;
}
