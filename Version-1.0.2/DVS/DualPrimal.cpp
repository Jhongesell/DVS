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


#include <stdlib.h>
#include "Definiciones.hpp"
#include "DualPrimal.hpp"




DualPrimal::DualPrimal(Interchange &inter)
{
   XP = NULL;
   YP = NULL;
   SP.resize(2);
   SP[0] = NULL;
   SP[1] = NULL;
   this->inter = &inter;
   //~ bdPrimals = inter.getBdPrimals();
   //~ bdDuals = inter.getBdDuals();
   //~ dualMult = inter.getDualMult();
   nD = inter.getND();
   nDual = inter.getNDuals();
   nP = inter.getNP();
   nPrimal = inter.getNPrimals();
   printf("\nprimals = %d %d   duals %d  %d", nP, nPrimal, nD, nDual);
   fflush(stdout);
   genMats();
}

void DualPrimal::a(ldouble *u, ldouble *v)
{
   /* v = a(u) is the average operator. u[nBd], v[nBd] are arrays of
      boundary values in DualP
   */
   int i, j;
   ldouble val, mult;
   for (i = 0; i < nD; i++)
   {
      mult = inter->bds->dualMult[i];
      val = 0.0;
      for (j = 0; j < mult; j++) val += u[inter->bds->bdDuals[i][j]->index];
      val /= mult;
      for (j = 0; j < mult; j++) v[inter->bds->bdDuals[i][j]->index] = val;
   }
}


void DualPrimal::j(ldouble *u, ldouble *v)
{
   /* v = j(u) is the jump operator. u[nBd], v[nBd] are arrays of
      boundary values in DualP
   */
   int i, j;
   ldouble val, mult;

   for (i = 0; i < nD; i++)
   {
      val = 0.0;
      mult =  inter->bds->dualMult[i];  // multiplicity
      for (j = 0; j < mult; j++) val += u[inter->bds->bdDuals[i][j]->index];
      val /= mult;
      for (j = 0; j < mult; j++) v[inter->bds->bdDuals[i][j]->index] = u[inter->bds->bdDuals[i][j]->index] - val;
   }
}

void DualPrimal::calcValues(ldouble *u)
{
   /* Computes all values of the solution using the dual values in u[]
      and puts the results into Subdomain scratch 0 array  */
   inter->clear(0);
   toSubdomains(0, u);
   inter->multOp(0, 1);
   inter->rhs(0);
   inter->diff(1, 0, 1);
   solveAPP(0, 1, 0, 2);
   toSubdomains(0, u);
   inter->knownValues(0);
}

void DualPrimal::fromSubdomains(int sc, ldouble *u)
{
   /* Transfers (dual) data from subdomains to global vector u[]   */
   inter->fromSubdomains(sc);
   size_t i, j;
   for (i = 0; i < inter->bds->bdDuals.size(); i++)
   {
      for (j = 0; j < inter->bds->bdDuals[i].size(); j++)
      {
         u[inter->bds->bdDuals[i][j]->index] = inter->rbdValues(inter->bds->bdDuals[i][j]->subd, inter->bds->bdDuals[i][j]->node);
      }
   }
}

void DualPrimal::genMats(void)
{
   ce.nameClassFunct("DualPrimal", "genMats");

   int i;
   //printf("\ngenMats %d",nPrimal);
   if (nPrimal > 0)
   {
      SP.resize(2);
      XP = new (nothrow) ldouble[nPrimal];
      if (XP == NULL) ce.memoryError("XP");
      YP = new (nothrow) ldouble[nPrimal];
      if (YP == NULL) ce.memoryError("YP");
      for (i = 0; i < nPrimal; i++) XP[0] = YP[i] = 0.0;
      SP[0] = inter->calcSP(0);
      SP[1] = inter->calcSP(1);
   }
}

void DualPrimal::multS(ldouble *u, ldouble *v)
{
   /* v = S*u  arrays u and v must be distinct */
   size_t i, j;
   int k;
   for (k = 0; k < nDual; k++) v[k] = 0.0;  // clear v
   inter->clear(0);
   toSubdomains(0, u);
   inter->multOp(0, 1);
   solveAPP(0, 1, 2, 3);
   inter->multOp(2, 0);
   inter->fromSubdomains(1);
   inter->diffValues(0);
   for (i = 0; i < inter->bds->bdDuals.size(); i++)
   {
      for (j = 0; j < inter->bds->bdDuals[i].size(); j++)
      {
         v[inter->bds->bdDuals[i][j]->index] += inter->rbdValues(inter->bds->bdDuals[i][j]->subd, inter->bds->bdDuals[i][j]->node);
      }
   }
}


void DualPrimal::solveAPP(int sp, int sc1, int sc2, int sc3)
{
   /* solve subdomains sc2 = APP-1(sc1). sc3 is an available scratch array */
   size_t p, i, j;
   int s;
   if (nPrimal == 0)
   {
      inter->inverse(sp, sc1, sc2);
      return;
   }
   for (s = 0; s < nPrimal; s++) YP[s] = 0.0;
   inter->inverse(sp, sc1, sc2);
   inter->multOp(sc2, sc2);
   inter->getPrimals(sc1);
   inter->diffValues(sc2);

   for (i = 0; i < inter->bds->bdPrimals.size(); i++)
   {
      p = inter->bds->bdPrimals[i][0]->index;
      for (j = 0; j < inter->bds->bdPrimals[i].size(); j++)
         YP[p] += inter->rbdValues(inter->bds->bdPrimals[i][j]->subd, inter->bds->bdPrimals[i][j]->node);
   }
   SP[sp]->solve(XP, YP);
   inter->clear(sc2);
   for (i = 0; i < inter->bds->bdPrimals.size(); i++)
   {
      p = inter->bds->bdPrimals[i][0]->index;
      for (j = 0; j < inter->bds->bdPrimals[i].size(); j++)
         inter->sbdValues(inter->bds->bdPrimals[i][j]->subd, inter->bds->bdPrimals[i][j]->node, XP[p]);
   }
   inter->setPrimals(sc2);
   inter->multOp(sc2, sc3);
   inter->diff(sc3, sc1, sc3);
   inter->inverse(sp, sc3, sc2);
   for (i = 0; i < inter->bds->bdPrimals.size(); i++)
   {
      p = inter->bds->bdPrimals[i][0]->index;
      for (j = 0; j < inter->bds->bdPrimals[i].size(); j++)
         inter->sbdValues(inter->bds->bdPrimals[i][j]->subd, inter->bds->bdPrimals[i][j]->node, XP[p]);
   }
   inter->setPrimals(sc2);
}

void DualPrimal::solveS(ldouble *u,  ldouble *v)
{
   /*    u = (S-1)(v)                        */
   size_t i, j, p;
   int k;
   inter->clear(0);
   toSubdomains(0, v);
   if (nPrimal == 0)
      inter->inverse(1, 0, 0);
   else
   {
      inter->inverse(1, 0, 1);
      inter->multOp(1, 1);
      for (k = 0; k < nPrimal; k++) YP[k] = 0.0;
      inter->getPrimals(1);
      for (i = 0; i < inter->bds->bdPrimals.size(); i++)
      {
         p = inter->bds->bdPrimals[i][0]->index;
         for (j = 0; j < inter->bds->bdPrimals[i].size(); j++) YP[p] -= inter->rbdValues(inter->bds->bdPrimals[i][j]->subd, inter->bds->bdPrimals[i][j]->node);
      }
      SP[1]->solve(XP, YP);
      inter->clear(1);
      for (i = 0; i < inter->bds->bdPrimals.size(); i++)
      {
         p = inter->bds->bdPrimals[i][0]->index;
         for (j = 0; j < inter->bds->bdPrimals[i].size(); j++) inter->sbdValues(inter->bds->bdPrimals[i][j]->subd, inter->bds->bdPrimals[i][j]->node, XP[p]);
      }
      inter->setPrimals(1);
      inter->multOp(1, 1);
      inter->diff(1, 0, 1);
      inter->inverse(1, 1, 0);
   }
   fromSubdomains(0, u);
}

void DualPrimal::toSubdomains(int sc,  ldouble *u)
{
   size_t i, j;
   for (i = 0; i < inter->bds->bdDuals.size(); i++)
      for (j = 0; j < inter->bds->bdDuals[i].size(); j++)
         inter->sbdValues(inter->bds->bdDuals[i][j]->subd, inter->bds->bdDuals[i][j]->node, u[inter->bds->bdDuals[i][j]->index]);
   inter->toSubdomains(sc);
}
