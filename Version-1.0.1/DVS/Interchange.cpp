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
#include "Interchange.hpp"
#include "BandSolve.hpp"
#include "BandCholesky.hpp"



Interchange::Interchange(PropDef &props)
{
   ce.nameClassFunct("Interchange", "Interchange");

   int i, j;

   this->props = &props;

   f = NULL;
   g = NULL;
   zero = NULL;
   one = NULL;

   genGeom();

   Interchange *inter = this;

   bds = new CreateBdNodes();
   if (bds == NULL) ce.memoryError("bds");
   initialize(nOmega);

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



void Interchange::initialize(int nOmega)
{
   ce.nameClassFunct("Interchange", "initialize");

   int i, j, k;
   size_t r, s;

   bds->maxBd = getMaxBdSize();
   bds->hbd.resize(bds->maxBd * nOmega + 1);

   vector<InternalBd*> xbds;
   for (i = 0; i < nOmega; i++)
   {
      xbds = omegas[i]->getInternalBd();
      for (r = 0; r < xbds.size(); r++) bds->hbd[++bds->ibd] = xbds[r];
   }

   bds->bdAll.resize((bds->maxBd * nOmega + 1));
   for (r = 0; r < bds->bdAll.size(); r++) bds->bdAll[r].resize(8);
   //~ printf("\n%d %d %d\n",maxBd,nOmega,bdAll.size());
   HeapSort hs(bds->hbd, bds->ibd);
   hs.sort();

   //~ printf("\nibd %d",ibd);


   int **pntype = new (nothrow) int*[nOmega];
   if (pntype == NULL) ce.memoryError("pntype");
   for (i = 0; i < nOmega; i++) pntype[i] = getNtype(i);
   i = 1;
   while (i <= bds->ibd)
   {
      for (j = 1; i + j <= bds->ibd; j++)
      {
         if (hs.rr(i)->equals(hs.rr(i + j))) continue;
         else break;
      }
      vector<BdNode*> bdSet(j);   // j = node multiplicity
      for (k = 0; k < j; k++)
      {
         bdSet[k] = new (nothrow) BdNode(hs.rr(i + k)->rsubd(), hs.rr(i + k)->rbd(), hs.rr(i + k)->rdp(), j);
         if (bdSet[k] == NULL) ce.memoryError("bdSet[k]", k);
         pntype[bdSet[k]->subd][hs.rr(i + k)->rnode()] |= j << 12;   // Assign multiplicity to ntype[subdomain][node]
      }
      bds->bdAll[bds->ibdAll++] = bdSet;
      i += j;
   }
   for (i = 0; i < nOmega; i++) setNtype(i, pntype[i]);
   delete []pntype;
   pntype = NULL;


   for (i = 0; i < bds->ibdAll; i++)
   {
      if (bds->bdAll[i][0]->index >= 0) bds->nD++;
      else bds->nP++;
   }
   //~ printf("\nibdAll %d\n",ibdAll);
   bds->bdDuals.resize(bds->nD);
   bds->bdPrimals.resize(bds->nP);
   bds->dualMult = new (nothrow) int[bds->nD];
   if (bds->dualMult == NULL) ce.memoryError("dualMult");
   for (i = 0; i < bds->nD; i++) bds->dualMult[i] = 0;
   bds->nP = bds->nD = 0;
   for (i = 0; i < bds->ibdAll; i++)
   {
      if (bds->bdAll[i][0]->index >= 0)
      {
         bds->dualMult[bds->nD] =  bds->bdAll[i][0]->mult;
         bds->bdDuals[bds->nD++] = bds->bdAll[i];
      }
      else bds->bdPrimals[bds->nP++] = bds->bdAll[i];
   }
   //~ printf("\nnD %d nP %d\n",nD, nP);
   bds->nDual = bds->nPrimal = 0;
   for (r = 0; r < bds->bdDuals.size(); r++)
      for (s = 0; s < bds->bdDuals[r].size(); s++) bds->bdDuals[r][s]->index = bds->nDual++;

   for (r = 0; r < bds->bdPrimals.size(); r++)
   {
      for (s = 0; s < bds->bdPrimals[r].size(); s++) bds->bdPrimals[r][s]->index = bds->nPrimal;
      bds->nPrimal++;
   }
}



// generate the goeometry
void  Interchange::genGeom(void)
{
   ce.nameClassFunct("Interchange", "genGeom");

   int i;
   zero = new (nothrow) Constant(0.0);
   if (zero == NULL) ce.memoryError("zero");
   one = new (nothrow) Constant(1.0);
   if (one == NULL) ce.memoryError("one");

   mesh = new (nothrow) int[6];
   if (mesh == NULL) ce.memoryError("mesh");
   mesh[0] = props->getInt("Nx", 3);
   mesh[1] = props->getInt("Ny");
   mesh[2] = props->getInt("Nz");
   mesh[3] = props->getInt("Mx");
   mesh[4] = props->getInt("My");
   mesh[5] = props->getInt("Mz");
   if (mesh[1] == 0) mesh[1] = mesh[0];
   if (mesh[2] == 0) mesh[2] = mesh[0];
   if (mesh[3] == 0) mesh[3] = mesh[0];
   if (mesh[4] == 0) mesh[4] = mesh[3];
   if (mesh[5] == 0) mesh[5] = mesh[3];
   nDim = props->getInt("d", 0);
   if (nDim < 2 || nDim > 3) ce.fatalError(1, "No se especifico la dimension del problema");

   swprint = props->getInt("p", 0);
   ldouble a[3];
   ldouble b[3];
   a[0] = props->getDouble("ax", 0.0);
   a[1] = props->getDouble("ay", 0.0);
   a[2] = props->getDouble("az", 0.0);
   b[0] = props->getDouble("bx", 0.0);
   b[1] = props->getDouble("by", 0.0);
   b[2] = props->getDouble("bz", 0.0);
   c = props->getDouble("c", 0.0);
   op = new (nothrow) EllipOp(nDim, a, b, c);
   if (op == NULL) ce.memoryError("op");

   domain = new (nothrow) ldouble *[3];
   if (domain == NULL) ce.memoryError("domain");
   for (i = 0; i < 3; i++)
   {
      domain[i] = new (nothrow) ldouble[2];
      if (domain[i] == NULL) ce.memoryError("domain[i]", i);
   }
   domain[0][0] = props->getDouble("Xmin", 0.0);
   domain[0][1] = props->getDouble("Xmax", 1.0);
   domain[1][0] = props->getDouble("Ymin", 0.0);
   domain[1][1] = props->getDouble("Ymax", 1.0);
   domain[2][0] = props->getDouble("Zmin", 0.0);
   domain[2][1] = props->getDouble("Zmax", 1.0);
   f = NULL;
   LookUpFunction lf;
   sf = props->getString("f", "");
   printf("\nPrint f %s", sf);
   if (sf[0] != 0) f = lf.getF(sf);
   fc = props->getDouble("fc", 0.0);
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
   sg = props->getString("g", "");
   if (sf[0] != 0) g = lf.getF(sg);
   gc = props->getDouble("gc", 0.0);
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
   op->setF(*f);
   op->setG(*g);

   prim = props->getString("pr", "noPrimal");
   if (strcmp(prim, "vertex") == 0) primal = new (nothrow) VertPrimal();
   else if (strcmp(prim, "vertedge") == 0) primal = new (nothrow) VertEdgePrimal();
   else if (strcmp(prim, "all") == 0) primal = new (nothrow) AllPrimal();
   else primal = new (nothrow) NoPrimal();
   if (primal == NULL) ce.memoryError("primal");

   method = props->getString("m", "MF1");
   printf("\nmethod %s", method);
   printf("\ndim= %d %d %d %d %d %d %d %s", nDim, mesh[0], mesh[1], mesh[2], mesh[3], mesh[4], mesh[5], method);
#ifdef __Double__
   printf("\nax = %f ay = %f bx = %f by = %f c = %f\n", a[0], a[1], b[0], b[1], c);
#else
   printf("\nax = %Lf ay = %Lf bx = %Lf by = %Lf c = %Lf\n", a[0], a[1], b[0], b[1], c);
#endif

   nOmega = (nDim == 1 ? mesh[0] : (nDim == 2 ? mesh[0] * mesh[1] : mesh[0] * mesh[1] * mesh[2]));
   omegas.resize(nOmega);
   for (i = 0; i < nOmega; i++)
   {
      omegas[i] = new (nothrow) RectSub(i, nDim, mesh, domain, *op, *primal);
      if (omegas[i] == NULL) ce.memoryError("omegas[i]", i);
   }
}

//~  Computes APP-1 where sp = 0 for APP with interior nodes
//~ sp = 1 for APP with interior + dual nodes
//~ This method is parallelizable to a point.
Solvable *Interchange::calcSP(int sp)
{
   ce.nameClassFunct("Interchange", "calcSP");

   int p, q, s;
   
   MatrizDispersa *SC;
   int ban = nDim == 2 ? 9 : bds->nPrimal/2;
   SC = new (nothrow) MatrizDispersa(bds->nPrimal, bds->nPrimal, ban, "SC");
   if (SC == NULL) ce.memoryError("SC");
   
   //~ ldouble **SC = new (nothrow) ldouble *[bds->nPrimal];
   //~ if (SC == NULL) ce.memoryError("SC");
   //~ for (p = 0; p < bds->nPrimal; p++)
   //~ {
      //~ SC[p] = new (nothrow) ldouble[bds->nPrimal];
      //~ if (SC[p] == NULL) ce.memoryError("SC[p]", p);
   //~ }
   //~ for (p = 0; p < bds->nPrimal; p++)
      //~ for (q = 0; q < bds->nPrimal; q++) SC[p][q] = 0.0;

   size_t i, j, k , r;
   ldouble xval;
   for (i = 0; i < bds->bdPrimals.size(); i++)
   {
      p = bds->bdPrimals[i][0]->index;
      for (j = 0; j < bds->bdPrimals[i].size(); j++)
      {
         s = bds->bdPrimals[i][j]->subd;
         calcula(s, bds->bdPrimals[i][j]->node, sp);
         for (k = 0; k < bds->bdPrimals.size(); k++)
         {
            q = bds->bdPrimals[k][0]->index;
            for (r = 0; r < bds->bdPrimals[k].size(); r++)
               if (bds->bdPrimals[k][r]->subd == s) {
                  xval = SC->retorna(p,q);
                  SC->asigna(p,q, xval + getValue(s, 1, 0, bds->bdPrimals[k][r]->node));
               }
         }
      }
   }
   
//~ SC->visualizaMatricesInternas();
   // Resolucion del sistema lineal por medio del metodo BanSolve
   Solvable *ss = new (nothrow) BandSolve(bds->nPrimal, SC);
   if (ss == NULL) ce.memoryError("ss in BandSolve");
   
   delete SC;
   SC = NULL;
   
   //~ BandSolve *ss = new (nothrow) BandSolve(bds->nPrimal, SC);
   //~ if (ss == NULL) ce.memoryError("ss");

   //~ for (int ii = 0; ii < bds->nPrimal; ii++) delete []SC[ii];
   //~ delete []SC;
   //~ SC = NULL;

   return ss;
}
