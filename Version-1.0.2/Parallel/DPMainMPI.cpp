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



#include "Definiciones.hpp"
#include "DPMainMPI.hpp"

/// Constructor de la clase
DPMainMPI::DPMainMPI(int id, int np, PropDef &props, EllipOp &op) : EsquemaMEMPI(id, np)
{
   ce.nameClassFunct("DPMainMPI", "DPMainMPI");

   int i;
   zero = new (nothrow) Constant(0.0);
   if (zero == NULL) ce.memoryError("zero");

   one = new (nothrow) Constant(1.0);
   if (one == NULL) ce.memoryError("one");

   mesh = new (nothrow) int[6];
   if (mesh == NULL) ce.memoryError("mesh");
   mesh[0] = props.getInt("Nx", 3);
   mesh[1] = props.getInt("Ny");
   mesh[2] = props.getInt("Nz");
   mesh[3] = props.getInt("Mx");
   mesh[4] = props.getInt("My");
   mesh[5] = props.getInt("Mz");
   if (mesh[1] == 0) mesh[1] = mesh[0];
   if (mesh[2] == 0) mesh[2] = mesh[0];
   if (mesh[3] == 0) mesh[3] = mesh[0];
   if (mesh[4] == 0) mesh[4] = mesh[3];
   if (mesh[5] == 0) mesh[5] = mesh[3];
   nDim = props.getInt("d", 0);
   if (nDim < 2 || nDim > 3) ce.fatalError(1, "No se especifico la dimension del problema");

   swprint = props.getInt("p", 0);
   ldouble a[3];
   ldouble b[3];
   a[0] = props.getDouble("ax", 0.0);
   a[1] = props.getDouble("ay", 0.0);
   a[2] = props.getDouble("az", 0.0);
   b[0] = props.getDouble("bx", 0.0);
   b[1] = props.getDouble("by", 0.0);
   b[2] = props.getDouble("bz", 0.0);
   c = props.getDouble("c", 0.0);
   //~ op = new EllipOp(nDim, a, b, c);
   this->op = &op;
   domain = new (nothrow) ldouble *[3];
   if (domain == NULL) ce.memoryError("domain");
   for (i = 0; i < 3; i++)
   {
      domain[i] = new (nothrow) ldouble[2];
      if (domain[i] == NULL) ce.memoryError("domain[i]", i);
   }
   domain[0][0] = props.getDouble("Xmin", 0.0);
   domain[0][1] = props.getDouble("Xmax", 1.0);
   domain[1][0] = props.getDouble("Ymin", 0.0);
   domain[1][1] = props.getDouble("Ymax", 1.0);
   domain[2][0] = props.getDouble("Zmin", 0.0);
   domain[2][1] = props.getDouble("Zmax", 1.0);
   f = NULL;
   LookUpFunction lf;
   sf = props.getString("f", "");
   printf("\nPrint f %s", sf);
   if (sf[0] != 0) f = lf.getF(sf);
   fc = props.getDouble("fc", 0.0);
   if (f != NULL) f->setVar(fc);
   else
   {
      f = new (nothrow) Constant(0.0);
      if (f == NULL) ce.memoryError("f");
   }
#ifdef __Double__
   printf("\nf  %f", fc);
   printf("\ngetVar %f", f->getVar());
#else
   printf("\nf  %Lf", fc);
   printf("\ngetVar %Lf", f->getVar());
#endif
   g = NULL;
   sg = props.getString("g", "");
   if (sf[0] != 0) g = lf.getF(sg);
   gc = props.getDouble("gc", 0.0);
#ifdef __Double__
   printf("\ng %f", gc);
#else
   printf("\ng %Lf", gc);
#endif
   if (g != NULL) g->setVar(gc);
   else
   {
      g = new (nothrow) Constant(0.0);
      if (g == NULL) ce.memoryError("g");
   }
#ifdef __Double__
   printf("\n%s %f %s %f", sf, fc, sg, gc);
#else
   printf("\n%s %Lf %s %Lf", sf, fc, sg, gc);
#endif
   op.setF(*f);
   op.setG(*g);
   prim = props.getString("pr", "noPrimal");
   if (strcmp(prim, "vertex") == 0) primal = new (nothrow) VertPrimal();
   else if (strcmp(prim, "vertedge") == 0) primal = new (nothrow) VertEdgePrimal();
   else if (strcmp(prim, "all") == 0) primal = new (nothrow) AllPrimal();
   else primal = new (nothrow) NoPrimal();
   if (primal == NULL) ce.memoryError("primal");

   method = props.getString("m", "CMG");
   printf("\nmethod %s", method);
   printf("\ndim= %d %d %d %d %d %d %d %s", nDim, mesh[0], mesh[1], mesh[2], mesh[3], mesh[4], mesh[5], method);
#ifdef __Double__
   printf("\nax = %f ay = %f bx = %f by = %f c = %f\n", a[0], a[1], b[0], b[1], c);
#else
   printf("\nax = %Lf ay = %Lf bx = %Lf by = %Lf c = %Lf\n", a[0], a[1], b[0], b[1], c);
#endif
   fflush(stdout);
   nOmega = (nDim == 1 ? mesh[0] : (nDim == 2 ? mesh[0] * mesh[1] : mesh[0] * mesh[1] * mesh[2]));


   // Controlador del esquema M-E
   if (id == 0)
   {
      time(&t1);
      printf("\nMaestro ID=%d\n", id);
      fflush(stdout);

      // Manda el numero de tareas que deben hacer los nodos esclavos
      generaRepartoCarga(nOmega);
   }
   else
   {
      // Numero de tareas a soportar
      MPI::COMM_WORLD.Recv(&nta, 1, MPI::INT, 0, 1);
      printf("\nEsclavo ID=%d  Tareas = %d\n", id, nta);
      fflush(stdout);
   }
}



/// Destructor de la clase
DPMainMPI::~DPMainMPI()
{
   int i;
   if (id == 0) // Maestro
   {
      time(&t2);
      printf("\nTiempo Calculo: %f\n", difftime(t2, t1));
      fflush(stdout);

   }
   else
   {
      for (i = 0; i < nta; i++) delete omegas[i];
   }

   delete zero;
   delete one;
   delete []mesh;
   for (i = 0; i < 3; i++) delete [] domain[i];
   delete []domain;


   delete f;
   delete g;
   delete []sf;
   delete []sg;
   delete []prim;
   delete []method;
   delete primal;
}


void DPMainMPI::deleteInternalBd(void)
{
   size_t i;
   for (i = 0; i < hbd.size(); i++)
   {
      InternalBd *abc = hbd[i];
      delete abc;
      abc = NULL;
   }
   hbd.clear();
}



/// Esclavo
void DPMainMPI::Esclavo(void)
{
   ce.nameClassFunct("DPMainMPI", "Esclavo");

   int indg, i, k, indl, sw;

   // Crea las geometrias correspondientes al nodo esclavo
   omegas.resize(nta);
   for (i = 0; i < nta; i++)
   {
      MPI::COMM_WORLD.Recv(&msa, 2, MPI::INT, 0, 1);
      indl = msa[0];
      indg = msa[1];
      //~ printf("\nLLego %d -> %d %d\n", i, indl, indg);
      //~ fflush(stdout);
      omegas[i] = new (nothrow) RectSub(indg, nDim, mesh, domain, *op, *primal);
      if (omegas[i] == NULL) ce.memoryError("omegas[i]", i);
   }

   sw = 1;
   while (sw)
   {

      MPI::COMM_WORLD.Recv(&msa, 10, MPI::INT, 0, 1);
      indl = msa[0];
      indg = msa[1];
      sw   = msa[2];
      //~ printf("\nLLego  %d %d (%d)\n", indg, indl,sw);
      //~ fflush(stdout);

      switch (sw)
      {
      case 0:
         printf("\nFinalizacion de tareas del esclavo %d\n", id);
         fflush(stdout);
         break;
      case 2:  // return omegas[e]->getNtype();
      {
         int tm = omegas[indl]->getNP();
         int *arr;
         mss[0] = tm;
         MPI::COMM_WORLD.Send(&mss, 1, MPI::INT, 0, 1);

         arr = omegas[indl]->getNtype();
         MPI::COMM_WORLD.Send(arr, tm, MPI::INT, 0, 1);
         break;
      }
      case 3:  // return omegas[e]->getInternalBd();
      {
         hbd = omegas[indl]->getInternalBd();
         int tm = hbd.size();
         mss[0] = tm;
         MPI::COMM_WORLD.Send(&mss, 1, MPI::INT, 0, 1);
         int s, n, b, i, d;
         ldouble *c = new (nothrow) ldouble[5];
         if (c == NULL) ce.memoryError("c");
         for (int j = 0; j < tm; j++)
         {
            hbd[j]->getval(s, n, b, i, d, c);
            mss[0] = s;
            mss[1] = n;
            mss[2] = b;
            mss[3] = i;
            mss[4] = d;
            MPI::COMM_WORLD.Send(&mss, 5, MPI::INT, 0, 1);
#ifdef __Double__
            MPI::COMM_WORLD.Send(c, d, MPI::DOUBLE, 0, 1);
#else
            MPI::COMM_WORLD.Send(c, d, MPI::LONG_DOUBLE, 0, 1);
#endif
         }
         delete []c;
         deleteInternalBd();
         break;
      }
      case 6:  // getValue(s,1,0,bdPrimals[k][r]->rnode());
      {
         ldouble val = (omegas[indl]->getValue(msa[3], msa[5]) - omegas[indl]->getValue(msa[4], msa[5]));
         //double val = omegas[indl]->getValue(msa[3], msa[4]);
#ifdef __Double__
         MPI::COMM_WORLD.Send(&val, 1, MPI::DOUBLE, 0, 1);
#else
         MPI::COMM_WORLD.Send(&val, 1, MPI::LONG_DOUBLE, 0, 1);
#endif
         break;
      }
      case 7:  // omegas[e]->diff(sc3, sc1, sc2);
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->diff(msa[3], msa[4], msa[5]);
         break;
      }
      case 9:  // omegas[e]->knownValues(sc);
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->knownValues(msa[3]);
         break;
      }
      case 11: // omegas[e]->rhs(sc);
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->rhs(msa[3]);
         break;
      }
      case 13: // omegas[e]->genInv(type);
      {
         omegas[indl]->genInv(msa[3]);
         break;
      }
      case 14: // omegas[e]->getCoordNode(n, x);
      {
         ldouble *arr = new (nothrow) ldouble[nDim];
         if (arr == NULL) ce.memoryError("arr");
         for (int j = 0; j < nDim; j++) arr[j] = 0.0;
         omegas[indl]->getCoordNode(msa[3], arr);
#ifdef __Double__
         MPI::COMM_WORLD.Send(arr, nDim, MPI::DOUBLE, 0, 1);
#else
         MPI::COMM_WORLD.Send(arr, nDim, MPI::LONG_DOUBLE, 0, 1);
#endif
         delete []arr;
         break;
      }
      case 15: // omegas[i]->print(s, sc);
      {
         char *cad = new (nothrow) char[msa[4]];
         if (cad == NULL) ce.memoryError("cad");
         MPI::COMM_WORLD.Recv(cad, msa[4], MPI::CHAR , 0, 1);
         int i;
         for (i = 0; i < nta; i++) omegas[i]->print(cad, msa[3]);
         fflush(stdout);
         delete []cad;
         break;
      }
      case 150: // omegas[i]->print(s, sc);
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->print(msa[3]);
         fflush(stdout);
         break;
      }
      case 16: // omegas[s]->diffValues(sc, bdValues[s]);
      {
         ldouble *arr = new (nothrow) ldouble[msa[4]];
         if (arr == NULL) ce.memoryError("arr");
         //for (int j = 0; j < msa[4]; j++) arr[j] = 0.0;
#ifdef __Double__
         MPI::COMM_WORLD.Recv(arr, msa[4], MPI::DOUBLE , 0, 1);
#else
         MPI::COMM_WORLD.Recv(arr, msa[4], MPI::LONG_DOUBLE , 0, 1);
#endif
         //~ for (int j = 0; j < msa[4]; j++) printf ("\nLlego E %f",arr[j]);
         //~ fflush(stdout);
         omegas[indl]->diffValues(msa[3], arr);
         //~ for (int j = 0; j < msa[4]; j++) printf ("\nEnvio E%f",arr[j]);
         //~ fflush(stdout);
#ifdef __Double__
         MPI::COMM_WORLD.Send(arr, msa[4], MPI::DOUBLE, 0, 1);
#else
         MPI::COMM_WORLD.Send(arr, msa[4], MPI::LONG_DOUBLE, 0, 1);
#endif
         delete []arr;
         break;
      }
      case 17: // omegas[s]->getValues(sc, bdValues[s]);
      {
         ldouble *arr = new (nothrow) ldouble[msa[4]];
         if (arr == NULL) ce.memoryError("arr");
         for (int j = 0; j < msa[4]; j++) arr[j] = 0.0;
         omegas[indl]->getValues(msa[3], arr);
         //~ for (int j = 0; j < msa[4]; j++) printf ("\nEnvio %f",arr[j]);
         //~ fflush(stdout);
#ifdef __Double__
         MPI::COMM_WORLD.Send(arr, msa[4], MPI::DOUBLE, 0, 1);
#else
         MPI::COMM_WORLD.Send(arr, msa[4], MPI::LONG_DOUBLE, 0, 1);
#endif
         delete []arr;
         break;
      }
      case 18: // omegas[s]->getPrimals(sc, bdValues[s]);
      {
         ldouble *arr = new (nothrow) ldouble[msa[4]];
         if (arr == NULL) ce.memoryError("arr");
         for (int j = 0; j < msa[4]; j++) arr[j] = 0.0;
         omegas[indl]->getPrimals(msa[3], arr);
         //~ for (int j = 0; j < msa[4]; j++) printf ("\nEnvio %f",arr[j]);
         //~ fflush(stdout);
#ifdef __Double__
         MPI::COMM_WORLD.Send(arr, msa[4], MPI::DOUBLE, 0, 1);
#else
         MPI::COMM_WORLD.Send(arr, msa[4], MPI::LONG_DOUBLE, 0, 1);
#endif
         delete []arr;
         break;
      }
      case 19: // omegas[s]->setPrimals(sc, bdValues[s]);
      {
         ldouble *arr = new (nothrow) ldouble[msa[4]];
         if (arr == NULL) ce.memoryError("arr");
         for (int j = 0; j < msa[4]; j++) arr[j] = 0.0;
#ifdef __Double__
         MPI::COMM_WORLD.Recv(arr, msa[4], MPI::DOUBLE , 0, 1);
#else
         MPI::COMM_WORLD.Recv(arr, msa[4], MPI::LONG_DOUBLE , 0, 1);
#endif
         omegas[indl]->setPrimals(msa[3], arr);
         delete []arr;
         break;
      }
      case 20: // omegas[s]->setValues(sc, bdValues[s]);
      {
         ldouble *arr = new (nothrow) ldouble[msa[4]];
         if (arr == NULL) ce.memoryError("arr");
         for (int j = 0; j < msa[4]; j++) arr[j] = 0.0;
#ifdef __Double__
         MPI::COMM_WORLD.Recv(arr, msa[4], MPI::DOUBLE , 0, 1);
#else
         MPI::COMM_WORLD.Recv(arr, msa[4], MPI::LONG_DOUBLE , 0, 1);
#endif
         //~ for (int j = 0; j < msa[4]; j++) printf ("\nLllego %f",arr[j]);
         //~ fflush(stdout);
         omegas[indl]->setValues(msa[3], arr);
         delete []arr;
         break;
      }
      case 50: // Se agrega para poder actualizar Ntype en paralelo (ya que no se pasan datos por referencia)
      {
         int tm = omegas[indl]->getNP();
         mss[0] = tm;
         int *arr = new (nothrow) int [tm];
         if (arr == NULL) ce.memoryError("arr");
         MPI::COMM_WORLD.Send(&mss, 1, MPI::INT, 0, 1);
         MPI::COMM_WORLD.Recv(arr, tm, MPI::INT, 0, 1);
         omegas[indl]->setNtype(arr);
         delete []arr;
         break;
      }
      case 51: // Se agrega para poder calcular el maximo de Bdsize en paralelo
      {
         int i, n = 0;
         mss[0] = 0;
         for (i = 0; i < nta; i++)
         {
            n = omegas[i]->getBdSize();
            if (n > mss[0]) mss[0] = n;
         }
         MPI::COMM_WORLD.Send(&mss, 1, MPI::INT, 0, 1);
         break;
      }
      case 52:  // omegas[e]->clear(sc); all subdomains
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->clear(msa[3]);
         break;
      }
      case 53: // omegas[e]->multOp(sc1, sc2);
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->multOp(msa[3], msa[4]);
         break;
      }
      case 54:  // omegas[e]->inverse(sp, sc1, sc2);
      {
         int i;
         for (i = 0; i < nta; i++) omegas[i]->inverse(msa[3], msa[4], msa[5]);
         break;
      }
      case 55:  // omegas[e]->inverse(sp, sc1, sc2);
      {
         omegas[indl]->clear(0);
         omegas[indl]->setValue(0, msa[3], 1.0);
         omegas[indl]->multOp(0, 1);
         omegas[indl]->inverse(msa[4], 1, 2);
         omegas[indl]->multOp(2, 0);
         break;
      }

      }
   }

}




