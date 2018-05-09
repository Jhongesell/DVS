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
#include "InterchangeMPI.hpp"
#include "BandSolve.hpp"
#include "CreateBdNodes.hpp"



InterchangeMPI::InterchangeMPI(PropDef &props, EsquemaMEMPI &me)
{
   ce.nameClassFunct("InterchangeMPI", "InterchangeMPI");

   f = NULL;
   g = NULL;
   zero = NULL;
   one = NULL;

   this->props = &props;
   ME = &me;

   genGeom();

   Interchange *inter = this;

   bds = new CreateBdNodes();
   if (bds == NULL) ce.memoryError("bds");
   initialize(nOmega);


   //~ nD = bds->getND();
   //~ nDual = bds->getNDuals();
   //~ nP = bds->getNP();
   //~ nPrimal = bds->getNPrimals();
   //~ maxBd = bds->getMaxBd();
   int i, j;

   bdValues = new (nothrow) ldouble *[nOmega];
   if (bdValues == NULL) ce.memoryError("bdValues");
   for (i = 0; i < nOmega; i++)
   {
      bdValues[i] = new (nothrow) ldouble[bds->maxBd];
      if (bdValues[i] == NULL) ce.memoryError("bdValues[i]", i);
   }
   for (i = 0; i < nOmega; i++)
      for (j = 0; j < bds->maxBd; j++) bdValues[i][j] = 0.0;
}


int InterchangeMPI::getMaxBdSize(void)
{
   int e, max = 0;
   //printf("\n 51 (%d) (%d) (%d)\n",0, 0, 0);
   //fflush(stdout);
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 51;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++)
   {
      MPI::COMM_WORLD.Recv(&msa, 1, MPI::INT, e, 1);
      if (msa[0] > max) max = msa[0];
   }
   return max;
   // return getMaxBdSize(); // get max Bdsize in all subdomains
}


int *InterchangeMPI::getNtype(int e)
{
   ce.nameClassFunct("InterchangeMPI", "getNtype");

   ME->reparteCargaTrabajo(xnp, indl, e);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 2;
   //~ printf("\n 26 (%d) (%d) (%d)\n",xnp, indl, e);
   //~ fflush(stdout);
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);

   MPI::COMM_WORLD.Recv(&msa, 1, MPI::INT, xnp, 1);
   int tm = msa[0];
   int *arr = new (nothrow) int[tm];
   if(arr == NULL) ce.memoryError("arr");
   MPI::COMM_WORLD.Recv(arr, tm, MPI::INT, xnp, 1);
   return arr;
   // return omegas[e]->getNtype();
}

// Se agrega para poder actualizar Ntype en paralelo (ya que no se pasan datos por referencia)
void InterchangeMPI::setNtype(int e, int * arr)
{
   ME->reparteCargaTrabajo(xnp, indl, e);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 50;
   //~ printf("\n 50 (%d) (%d) (%d)\n",xnp, indl, e);
   //~ fflush(stdout);
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);

   MPI::COMM_WORLD.Recv(&msa, 1, MPI::INT, xnp, 1);
   int tm = msa[0];
   MPI::COMM_WORLD.Send(arr, tm, MPI::INT, xnp, 1);
}

vector<InternalBd*> InterchangeMPI::getInternalBd(int e)
{
   ce.nameClassFunct("InterchangeMPI", "getInternalBd");

   ME->reparteCargaTrabajo(xnp, indl, e);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 3;
   //~ printf("\n 27 (%d) (%d) (%d)\n",xnp, indl, e);
   //~ fflush(stdout);
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);

   MPI::COMM_WORLD.Recv(&msa, 1, MPI::INT, xnp, 1);
   int m, j, np = msa[0];
   vector<InternalBd*> *intBd = new vector<InternalBd*>(np);
   vector<InternalBd*> &pintBd = *intBd;
   int s, n, b, i, d;
   for (m = 0; m < np; m++)
   {
      MPI::COMM_WORLD.Recv(&msa, 5, MPI::INT, xnp, 1);
      s = msa[0];
      n = msa[1];
      b = msa[2];
      i = msa[3];
      d = msa[4];
      ldouble *xc = new (nothrow) ldouble[d];
      if (xc == NULL) ce.memoryError("xc");
#ifdef __Double__
      MPI::COMM_WORLD.Recv(xc, d, MPI::DOUBLE, xnp, 1);
#else
      MPI::COMM_WORLD.Recv(xc, d, MPI::LONG_DOUBLE, xnp, 1);
#endif
      pintBd[m] = new (nothrow) InternalBd(s, n, b, i, d, xc);
      if (pintBd[m] == NULL) ce.memoryError("pintBd[m]", m);
   }
   return pintBd;
   // return omegas[e]->getInternalBd();
}

void InterchangeMPI::calcula(int e, int node, int sp)
{
   ME->reparteCargaTrabajo(xnp, indl, e);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 55;
   mss[3] = node;
   mss[4] = sp;
   //~ printf("\n 55 (%d) (%d) (%d)\n",0, 0, 0);
   //~ fflush(stdout);
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
   //~ clear(e,0);
   //~ setValue(e, 0, node, 1.0);
   //~ multOp(e, 0, 1);
   //~ inverse(e, sp, 1, 2);
   //~ multOp(e, 2, 0);
}



// Clear scr[sc][] in all subdomains
void InterchangeMPI::clear(int sc)
{
   int e;
   //~ printf("\n 52 (%d) (%d) (%d)\n",0, 0, 0);
   //~ fflush(stdout);
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 52;
   mss[3] = sc;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->clear(sc);
}


ldouble InterchangeMPI::getValue(int e, int scr1, int scr2, int node)
{
   ME->reparteCargaTrabajo(xnp, indl, e);
   //~ printf("\n 4 (%d) (%d) (%d)\n",xnp, indl, e);
   //~ fflush(stdout);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 6;
   mss[3] = scr1;
   mss[4] = scr2;
   mss[5] = node;
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
   ldouble val;
#ifdef __Double__
   MPI::COMM_WORLD.Recv(&val, 1, MPI::DOUBLE, xnp, 1);
#else
   MPI::COMM_WORLD.Recv(&val, 1, MPI::LONG_DOUBLE, xnp, 1);
#endif
   return val;
   // return getValue(s,1,0,bdPrimals[k][r]->rnode());
   // Que remplaza a estas otras
   // getValue(s, 1, bdPrimals[k][r]->rnode()) - getValue(s, 0, bdPrimals[k][r]->rnode());
}

// scr[sc3][] = scr[sc1][] - scr[sc2][] in all subdomains
void InterchangeMPI::diff(int sc3, int sc1, int sc2)
{
   //~ printf("\n 5 (%d) (%d) (%d)\n",0, 0, 0);
   //~ fflush(stdout);
   int e;
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 7;
   mss[3] = sc3;
   mss[4] = sc1;
   mss[5] = sc2;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->diff(sc3, sc1, sc2);
}

// scr[sc2][] = A(sp)-1(scr[sc1][])
void InterchangeMPI::inverse(int sp, int sc1, int sc2)
{
   int e;
   //~ printf("\n 6 (%d) (%d) (%d)\n",0, 0 0);
   //~ fflush(stdout);
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 54;
   mss[3] = sp;
   mss[4] = sc1;
   mss[5] = sc2;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->inverse(sp, sc1, sc2);
}


// scr[sc][] = Dirichlet boundary values of all subdomains
void InterchangeMPI::knownValues(int sc)
{
   //~ printf("\n 8 (%d) (%d) (%d)\n",0,0,0);
   //~ fflush(stdout);
   int e;
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 9;
   mss[3] = sc;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->knownValues(sc);
}

// scr[s2][] = A(scr[sc1][])
void InterchangeMPI::multOp(int sc1, int sc2)
{
   //~ printf("\n 9 (%d) (%d) (%d)\n",0,0,0);
   //~ fflush(stdout);
   int e;
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 53;
   mss[3] = sc1;
   mss[4] = sc2;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->multOp(sc1, sc2);
}

// scr[sc][] = initial right-hand-side (all subdomains)
void InterchangeMPI::rhs(int sc)
{
   //~ printf("\n 11 (%d) (%d) (%d)\n",0,0,0);
   //~ fflush(stdout);
   int e;
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 11;
   mss[3] = sc;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++) MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->rhs(sc);
}


void InterchangeMPI::genInv(int e, int type)
{
   ME->reparteCargaTrabajo(xnp, indl, e);
   //~ printf("\n 13 (%d) (%d) (%d)\n",xnp, indl, e);
   //~ fflush(stdout);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 13;
   mss[3] = type;
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
   // omegas[e]->genInv(type);
}

void InterchangeMPI::getCoordNode(int e, int n, ldouble *x)
{
   ME->reparteCargaTrabajo(xnp, indl, e);
   //~ printf("\n 14 (%d) (%d) (%d)\n",xnp, indl, e);
   //~ fflush(stdout);
   mss[0] = indl;
   mss[1] = e;
   mss[2] = 14;
   mss[3] = n;
   MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
#ifdef __Double__
   MPI::COMM_WORLD.Recv(x, nDim, MPI::DOUBLE, xnp, 1);
#else
   MPI::COMM_WORLD.Recv(x, nDim, MPI::LONG_DOUBLE, xnp, 1);
#endif
   // omegas[e]->getCoordNode(n, x);
}



void InterchangeMPI::print(const char *s, int sc)
{
   int e;
   //~ printf("\n 15 (%d) (%d) (%d)\n",0, 0, 0);
   //~ fflush(stdout);
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 15;
   mss[3] = sc;
   mss[4] = strlen(s) + 1;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++)
   {
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
      MPI::COMM_WORLD.Send(s, mss[4], MPI::CHAR, e, 1);
   }
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->print(s, sc);
}


void InterchangeMPI::print(int sc)
{
   int e;
   //~ printf("\n 15 (%d) (%d) (%d)\n",0, 0, 0);
   //~ fflush(stdout);
   mss[0] = 0;
   mss[1] = 0;
   mss[2] = 150;
   mss[3] = sc;
   //~ mss[4] = strlen(s)+1;
   for (e = 1; e <= ME->numeroProcesadoresUsar(); e++)
   {
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, e, 1);
      //~ MPI::COMM_WORLD.Send(s, mss[4], MPI::CHAR, e, 1);
   }
   // for (size_t i = 0; i < nOmega; i++) omegas[i]->print(sc);
}

/// bdValues[][] -= scr[sc][] in all subdomains
void InterchangeMPI::diffValues(int sc)
{
   for (int e = 0; e < nOmega; e++)
   {
      ME->reparteCargaTrabajo(xnp, indl, e);
      //~ printf("\n 17 (%d) (%d) (%d)\n",xnp, indl, e);
      //~ fflush(stdout);
      mss[0] = indl;
      mss[1] = e;
      mss[2] = 16;
      mss[3] = sc;
      mss[4] = bds->maxBd;
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
      //~ for (int j = 0; j < maxBd; j++) printf ("\nEnvio M%f",bdValues[e][j]);
      //~ fflush(stdout);
#ifdef __Double__
      MPI::COMM_WORLD.Send(bdValues[e], bds->maxBd, MPI::DOUBLE, xnp, 1);
      MPI::COMM_WORLD.Recv(bdValues[e], bds->maxBd, MPI::DOUBLE, xnp, 1);
#else
      MPI::COMM_WORLD.Send(bdValues[e], bds->maxBd, MPI::LONG_DOUBLE, xnp, 1);
      MPI::COMM_WORLD.Recv(bdValues[e], bds->maxBd, MPI::LONG_DOUBLE, xnp, 1);
#endif
      //~ for (int j = 0; j < maxBd; j++) printf ("\nLlego M%f",bdValues[e][j]);
      //~ fflush(stdout);
   }
   // for (int s = 0; s < nOmega; s++) omegas[s]->diffValues(sc, bdValues[s]);
}

// bdValues[][] = scr[sc][] from all subdomains
void InterchangeMPI::fromSubdomains(int sc)
{
   for (int e = 0; e < nOmega; e++)
   {
      ME->reparteCargaTrabajo(xnp, indl, e);
      //~ printf("\n 18 (%d) (%d) (%d)\n",xnp, indl, e);
      //~ fflush(stdout);
      mss[0] = indl;
      mss[1] = e;
      mss[2] = 17;
      mss[3] = sc;
      mss[4] = bds->maxBd;
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
#ifdef __Double__
      MPI::COMM_WORLD.Recv(bdValues[e], bds->maxBd, MPI::DOUBLE, xnp, 1);
#else
      MPI::COMM_WORLD.Recv(bdValues[e], bds->maxBd, MPI::LONG_DOUBLE, xnp, 1);
#endif
      //~ for (int j = 0; j < maxBd; j++) printf ("\nLlego %f",bdValues[e][j]);
      //~ fflush(stdout);
   }
   //for (int s = 0; s < nOmega; s++) omegas[s]->getValues(sc, bdValues[s]);
}

// bdValues[][] (primals only) = scr[sc][] (primals)
void InterchangeMPI::getPrimals(int sc)
{
   for (int e = 0; e < nOmega; e++)
   {
      ME->reparteCargaTrabajo(xnp, indl, e);
      //~ printf("\n 19 (%d) (%d) (%d)\n",xnp, indl, e);
      //~ fflush(stdout);
      mss[0] = indl;
      mss[1] = e;
      mss[2] = 18;
      mss[3] = sc;
      mss[4] = bds->maxBd;
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
#ifdef __Double__
      MPI::COMM_WORLD.Recv(bdValues[e], bds->maxBd, MPI::DOUBLE, xnp, 1);
#else
      MPI::COMM_WORLD.Recv(bdValues[e], bds->maxBd, MPI::LONG_DOUBLE, xnp, 1);
#endif
      //~ for (int j = 0; j < maxBd; j++) printf ("\nLlego %f",bdValues[e][j]);
      //~ fflush(stdout);
   }
   // for (int s = 0; s < nOmega; s++) omegas[s]->getPrimals(sc, bdValues[s]);
}

// scr[sc][] = bdValues  all subdomains
void InterchangeMPI::setPrimals(int sc)
{
   for (int e = 0; e < nOmega; e++)
   {
      ME->reparteCargaTrabajo(xnp, indl, e);
      //~ printf("\n 20 (%d) (%d) (%d)\n",xnp, indl, e);
      //~ fflush(stdout);
      mss[0] = indl;
      mss[1] = e;
      mss[2] = 19;
      mss[3] = sc;
      mss[4] = bds->maxBd;
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
#ifdef __Double__
      MPI::COMM_WORLD.Send(bdValues[e], bds->maxBd, MPI::DOUBLE, xnp, 1);
#else
      MPI::COMM_WORLD.Send(bdValues[e], bds->maxBd, MPI::LONG_DOUBLE, xnp, 1);
#endif
   }
   // for (int s = 0; s < nOmega; s++) omegas[s]->setPrimals(sc, bdValues[s]);
}

// scr[sc][] = bdValues[][]  all subdomains
void InterchangeMPI::toSubdomains(int sc)
{
   for (int e = 0; e < nOmega; e++)
   {
      ME->reparteCargaTrabajo(xnp, indl, e);
      //~ printf("\n 21 (%d) (%d) (%d)\n",xnp, indl, e);
      //~ fflush(stdout);
      mss[0] = indl;
      mss[1] = e;
      mss[2] = 20;
      mss[3] = sc;
      mss[4] = bds->maxBd;
      MPI::COMM_WORLD.Send(&mss, 10, MPI::INT, xnp, 1);
      //~ for (int j = 0; j < maxBd; j++) printf ("\nenvio %f",bdValues[e][j]);
      //~ fflush(stdout);
#ifdef __Double__
      MPI::COMM_WORLD.Send(bdValues[e], bds->maxBd, MPI::DOUBLE, xnp, 1);
#else
      MPI::COMM_WORLD.Send(bdValues[e], bds->maxBd, MPI::LONG_DOUBLE, xnp, 1);
#endif
   }
   // for (int s = 0; s < nOmega; s++) omegas[s]->setValues(sc, bdValues[s]);
}
