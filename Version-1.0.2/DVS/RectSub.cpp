// Grupo de Modelacion Matematica y Computacional
// Instituto de Goefisica
// Universidad Nacional Autonoma de Mexico
// Mexico D.F.
//
// Autores:
// Ismael Herrera Revilla
// Robert Yates Smit
// Antonio Carrillo Ledesma
// Alberto Rosas Medina
//
// Codigo liberado bajo la licencia GPL ver. 2.0



#include <stdlib.h>
#include "Definiciones.hpp"
#include "RectSub.hpp"
#include "EllipOp.hpp"
#ifndef GMM
#include "ICGM.hpp"
#include "IDQGMRES.hpp"
#include "BandSolve.hpp"
#include "BandCholesky.hpp"
#endif



RectSub::RectSub(int id, int nDim, int *mesh, ldouble **dom, EllipOp &op, Primal &primal)
{
   ce.nameClassFunct("RectSub", "RectSub");
   int i, j;


   this->op = &op;
   this->id = id;
   this->nDim = nDim;
   this->mesh = mesh;
   domain = new (nothrow) ldouble *[nDim];
   if (domain == NULL) ce.memoryError("domain");
   for (i = 0; i < nDim; i++)
   {
      domain[i] = new (nothrow) ldouble[2];
      if (domain[i] == NULL) ce.memoryError("domain[i]", i);
   }
   for (i = 0; i < nDim; i++)
      for (j = 0; j < 2; j++) domain[i][j] = 0.0;

   coord = new (nothrow) int[3];
   if (coord == NULL) ce.memoryError("coord");
   for (i = 0; i < 3; i++) coord[i] = 0;
   coordN = new (nothrow) int[3];
   if (coordN == NULL) ce.memoryError("coordN");
   for (i = 0; i < 3; i++) coordN[i] = 0;
   x = new (nothrow) ldouble[nDim];
   if (x == NULL) ce.memoryError("x");
   for (i = 0; i < nDim; i++) x[i] = 0.0;
   N = new (nothrow) int[3];
   if (N == NULL) ce.memoryError("N");
   for (i = 0; i < 3; i++) N[i] = 0;
   M = new (nothrow) int[3];
   if (M == NULL) ce.memoryError("M");
   for (i = 0; i < 3; i++) M[i] = 0;
   M1 = new (nothrow) int[3];
   if (M1 == NULL) ce.memoryError("M1");
   for (i = 0; i < 3; i++) M1[i] = 0;
   h = new (nothrow) ldouble[nDim];
   if (h == NULL) ce.memoryError("h");
   for (i = 0; i < nDim; i++) h[i] = 0.0;

   np = 1;
   hfac = 1.0;
   for (i = 0; i < nDim; i++)
   {
      N[i] = (i == 0 ? mesh[i] : mesh[i - 1] * N[i - 1]);
      M[i] = mesh[i + 3];
      M1[i] = (i == 0 ? M[i] + 1 : M1[i - 1] * (M[i] + 1));
      np *= (M[i] + 1);
      //printf("\n i(%d)  N=%d M=%d M1=%d np=%d",i, N[i],M[i],M1[i],np);
   }
   genNcoord(id, coordN, N);
   ldouble d;
   for (i = 0; i < nDim; i++)
   {
      d = (dom[i][1] - dom[i][0]) / mesh[i];
      domain[i][0] = dom[i][0] + coordN[i] * d;
      domain[i][1] = domain[i][0] + d;
      h[i] = (domain[i][1] - domain[i][0]) / M[i];
      hfac *= h[i];
      //printf("\n%d %d %1.15lf %1.15lf %1.15lf",id, i, domain[i][0], domain[i][1], h[i]);
   }
   bFloat = true;
   nBd = 0;
   ntype = new (nothrow) int[np];
   if (ntype == NULL) ce.memoryError("ntype");
   for (i = 0; i < np; i++) ntype[i] = 0;

   genNtype(primal);
   //for (i = 0; i < np; i++) printf("\n %d  %d",i,ntype[i]);

   for (i = 0; i < np; i++)
   {
      if ((ntype[i] & KNOWN) != 0 || (ntype[i] & PRIMAL) != 0) bFloat = false;
      if ((ntype[i] & INTBD) != 0) nBd++;
   }
   bdMap = new (nothrow) int[nBd];
   if (bdMap == NULL) ce.memoryError("bdMap");
   for (i = 0; i < nBd; i++) bdMap[i] = 0;

   nBd = 0;
   for (i = 0; i < np; i++) if ((ntype[i] & INTBD) != 0) bdMap[nBd++] = i;
   //printf("\nid %d",id);
   //for (i = 0; i < nBd; i++) printf("\n %d  %d",i,bdMap[i]);

   scr = new (nothrow) ldouble *[4];
   if (scr == NULL) ce.memoryError("scr");
   for (i = 0; i < 4; i++)
   {
      scr[i] = new (nothrow) ldouble[np];
      if (scr[i] == NULL) ce.memoryError("scr[i]", i);
   }
   for (i = 0; i < 4; i++)
      for (j = 0; j < np; j++) scr[i][j] = 0.0;

   coef = new (nothrow) ldouble[7];
   if (coef == NULL) ce.memoryError("coef");
   for (i = 0; i < 7; i++) coef[i] = 0.0;

#ifndef GMM
   X = new (nothrow) ldouble[np];
   if (X == NULL) ce.memoryError("X");
   for (i = 0; i < np; i++) X[i] = 0.0;

   Y = new (nothrow) ldouble[np];
   if (Y == NULL) ce.memoryError("Y");
   for (i = 0; i < np; i++) Y[i] = 0.0;
#endif
   genCoef(op);
}


void RectSub::clear(int s)
{
   for (int i = 0; i < np; i++) scr[s][i] = 0.0;
}

void RectSub::diff(int sc3, int sc1, int sc2)
{
   for (int i = 0; i < np; i++) scr[sc3][i] = scr[sc1][i] - scr[sc2][i];
}

// Si el problema es de coeficientes constantes aqui se genera el estencil
void RectSub::genCoef(EllipOp &op)
{
   ce.nameClassFunct("RectSub", "genCoef");
   int i;

#ifndef ESTABILIZA
   ldouble val;
#else
   ldouble val;
   ldouble PeX, PeY, PeZ, ExpX, ExpY, ExpZ;
#endif

   ldouble *a = op.getA();
   if (a == NULL) ce.memoryError("a");
   for (i = 0; i < 3; i++) Ca[i] = a[i];

   ldouble *b = op.getB();
   if (b == NULL)
   {
      b = new (nothrow) ldouble[3];
      if (b == NULL) ce.memoryError("b");
      for (i = 0; i < 3; i++) b[i] = 0.0;
   }
   for (i = 0; i < 3; i++) Cb[i] = b[i];


   ldouble c = op.getC();
   Cc = c;
   bsym = op.isSymmetric();

// Se asume coeficientes constantes
#ifdef COEFICIENTES_CONSTANTES

#ifndef ESTABILIZA
   // Se asume que no hay estabilizcion
   val = Ca[0] / (h[0] * h[0]);
   if (nDim > 1) val += Ca[1] / (h[1] * h[1]);
   if (nDim > 2) val += Ca[2] / (h[2] * h[2]);
   coef[0] = Cc + 2.0 * val; //coeficientes Uij
   coef[1] = (1.0 / h[0]) * (Cb[0] / 2.0 - Ca[0] / h[0]); //coeficientes ui+1,j
   coef[2] = -(1.0 / h[0]) * (Cb[0] / 2.0 + Ca[0] / h[0]); // Coeficientes ui-1,j
   if (nDim > 1)
   {
      coef[3] = (1.0 / h[1]) * (Cb[1] / 2.0 - Ca[1] / h[1]); // Coeficientes ui,j+1
      coef[4] = -(1.0 / h[1]) * (Cb[1] / 2.0 + Ca[1] / h[1]); // Coeficientes ui,j-1
   }
   if (nDim > 2)
   {
      coef[5] = (1.0 / h[2]) * (Cb[2] / 2.0 - Ca[2] / h[2]); // coeficientes ui,j,k+1
      coef[6] = -(1.0 / h[2]) * (Cb[2] / 2.0 + Ca[2] / h[2]); // coeficientes ui,j,k-1
   }
#else

   for (i = 0; i < 7; i++) coef[i]  = 0.0;

   val = Ca[0]/(h[0] * h[0]); // a0/hx^2
   if (nDim > 1) val += Ca[1]/(h[1] * h[1]); //  a0/(hx^2)+a1/(hy^2)
   if (nDim > 2) val += Ca[2]/(h[2] * h[2]); // a0/(hx^2)+a1/(hy^2)+a2/(hz^2)
   coef[0] = Cc + 2.0 * val + fabs(Cb[0])/h[0] + fabs(Cb[1])/h[1] //coeficientes Uij
   if (nDim > 2) coef[0] += fabs(Cb[2])/h[2]; // Para 3 dimensiones
   coef[1] = Cb[0]/(2.0*h[0]) - Ca[0]/(h[0] * h[0]) -fabs(Cb[0])/(2.0*h[0]); //coeficientes ui+1,j
   coef[2] = -Cb[0]/(2.0*h[0]) - Ca[0]/(h[0] * h[0]) -fabs(Cb[0])/(2.0*h[0]) ; // Coeficientes ui-1,j
   if (nDim > 1)
   {
      coef[3] = Cb[1]/(2.0*h[1]) - Ca[1]/(h[1] * h[1]) -fabs(Cb[1])/(2.0*h[1]); // Coeficientes ui,j+1
      coef[4] = -Cb[1]/(2.0*h[1]) - Ca[1]/(h[1] * h[1]) -fabs(Cb[1])/(2.0*h[1]); // Coeficientes ui,j-1
   }
   if (nDim > 2)
   {
      coef[5] = Cb[2]/(2.0*h[2]) - Ca[2]/(h[2] * h[2]) -fabs(Cb[2])/(2.0*h[2]); // coeficientes ui,j,k+1
      coef[6] = -Cb[2]/(2.0*h[2]) - Ca[2]/(h[2] * h[2]) -fabs(Cb[2])/(2.0*h[2]); // coeficientes ui,j,k-1
   }   
  
#endif

   //~ for (i = 0; i < 7; i++) printf("  %f ",coef[i]);
     //~ exit(1); 


   for (i = 0; i < 7; i++) coef[i] *= hfac;

#endif
}


// Si el problema es de coeficientes variables aqui se genera el estencil
void RectSub::genCoefVar(int ren)
{
   int i;
   ldouble x[3];
   ldouble val;
   ldouble PeX, PeY, PeZ, ExpX, ExpY, ExpZ;


   // Obtiene las coordenadas del nodo
   getCoord(ren, x);


   //~ // Coeficientes variables en A
   //~ if (x[0] < 0.0 || x[1] < 0.0) Ca[0] = 0.3, Ca[1] = 0.3;
   //~ else Ca[0] = 1.0, Ca[1] = 1.0;

#ifdef NoSimetricoVariable
   Cb[0] = (1.0 + x[1]) / 2.0;
   Cb[1] = 0.0;
   Cb[2] = 0.0;
#endif

#ifdef NoSimetricoNeumann
   Cb[0] = 1.0 / 2.0 * ( (1.0 - x[0] * x[0]) * (1.0 + x[1]));
   Cb[1] = -1.0 / 2.0 * (4.0 - ((1.0 + x[1]) * (1.0 + x[1])));
   Cb[2] = 0.0;
#endif

#ifdef NoSimetricoRotacional
   Cb[0] = x[1];
   Cb[1] = -x[0];
   Cb[2] = 0.0;
#endif


#ifndef ESTABILIZA
   val = Ca[0] / (h[0] * h[0]);
   if (nDim > 1) val += Ca[1] / (h[1] * h[1]);
   if (nDim > 2) val += Ca[2] / (h[2] * h[2]);
   coef[0] = Cc + 2.0 * val; //coeficientes Uij
   coef[1] = (1.0 / h[0]) * (Cb[0] / 2.0 - Ca[0] / h[0]); //coeficientes ui+1,j
   coef[2] = -(1.0 / h[0]) * (Cb[0] / 2.0 + Ca[0] / h[0]); // Coeficientes ui-1,j
   if (nDim > 1)
   {
      coef[3] = (1.0 / h[1]) * (Cb[1] / 2.0 - Ca[1] / h[1]); // Coeficientes ui,j+1
      coef[4] = -(1.0 / h[1]) * (Cb[1] / 2.0 + Ca[1] / h[1]); // Coeficientes ui,j-1
   }
   if (nDim > 2)
   {
      coef[5] = (1.0 / h[2]) * (Cb[2] / 2.0 - Ca[2] / h[2]); // coeficientes ui,j,k+1
      coef[6] = -(1.0 / h[2]) * (Cb[2] / 2.0 + Ca[2] / h[2]); // coeficientes ui,j,k-1
   }
#else
   for (i = 0; i < 7; i++) coef[i]  = 0.0;

   val = Ca[0]/(h[0] * h[0]); // a0/hx^2
   if (nDim > 1) val += Ca[1]/(h[1] * h[1]); //  a0/(hx^2)+a1/(hy^2)
   if (nDim > 2) val += Ca[2]/(h[2] * h[2]); // a0/(hx^2)+a1/(hy^2)+a2/(hz^2)
   coef[0] = Cc + 2.0 * val + fabs(Cb[0])/h[0] + fabs(Cb[1])/h[1]; //coeficientes Uij
   if (nDim > 2) coef[0] += fabs(Cb[2])/h[2]; // Para 3 dimensiones
   coef[1] = Cb[0]/(2.0*h[0]) - Ca[0]/(h[0] * h[0]) -fabs(Cb[0])/(2.0*h[0]); //coeficientes ui+1,j
   coef[2] = -Cb[0]/(2.0*h[0]) - Ca[0]/(h[0] * h[0]) -fabs(Cb[0])/(2.0*h[0]) ; // Coeficientes ui-1,j
   if (nDim > 1)
   {
      coef[3] = Cb[1]/(2.0*h[1]) - Ca[1]/(h[1] * h[1]) -fabs(Cb[1])/(2.0*h[1]); // Coeficientes ui,j+1
      coef[4] = -Cb[1]/(2.0*h[1]) - Ca[1]/(h[1] * h[1]) -fabs(Cb[1])/(2.0*h[1]); // Coeficientes ui,j-1
   }
   if (nDim > 2)
   {
      coef[5] = Cb[2]/(2.0*h[2]) - Ca[2]/(h[2] * h[2]) -fabs(Cb[2])/(2.0*h[2]); // coeficientes ui,j,k+1
      coef[6] = -Cb[2]/(2.0*h[2]) - Ca[2]/(h[2] * h[2]) -fabs(Cb[2])/(2.0*h[2]); // coeficientes ui,j,k-1
   } 
#endif

   //~ for (i = 0; i < 7; i++) printf("  %f ",coef[i]);
     //~ exit(1);

   for (i = 0; i < 7; i++) coef[i] *= hfac;
}


void RectSub::genInv(int type)
{
   ce.nameClassFunct("RectSub", "genInv");
   int i, nDual = 0;
   ldouble fac = 0.0;

   mapInt = new (nothrow) int[np];
   if (mapInt == NULL) ce.memoryError("mapInt");
   for (i = 0; i < np; i++) mapInt[i] = 0;

   nInt = 0;
   for (i = 0; i < np; i++) mapInt[i] = (isInterior(i) ? nInt++ : -1);
   mapFull = new (nothrow) int[np];
   if (mapFull == NULL) ce.memoryError("mapFull");
   for (i = 0; i < np; i++) mapFull[i] = 0;

   nFull = 0;
   for (i = 0; i < np; i++)
   {
      mapFull[i] = (isKnown(i) || isPrimal(i) ? -1 : nFull++);
      if (isDual(i)) nDual++;
   }

#ifdef GMM
   genInverse(mapInt, 0.0, 0);
   fac = ((type == 0 || !bFloat) ? 0.0 : 1.0 / nDual);
   genInverse(mapFull, fac, 1);
   
   //~ std::cout << "Matriz MapInt"<< MapInt << endl;
   //~ std::cout << "Matriz MapFull"<< MapFull << endl;
#else   
   inv.resize(2);
   inv[0] = genInverse(mapInt, 0.0);
   fac = ((type == 0 || !bFloat) ? 0.0 : 1.0 / nDual);
   inv[1] = genInverse(mapFull, fac);
#endif
}


#ifdef GMM
void RectSub::genInverse(int *map, ldouble fac, int ind)
{
   ce.nameClassFunct("RectSub", "genInverse");

   int b2 = 0, b1 = 0, Mx = M[0] + 1, Mxy = 0, i = 0, j = 0, k = 0, n = 0, m = 0, m1 = 0;
   Mxy = (nDim > 1 ? (M[0] + 1) * (M[1] + 1) : 0);
   for (int r = 0; r < np; r++) if (map[r] >= 0) n++;

   

   gmm::row_matrix< gmm::rsvector<double> > A(n, n);
   
   for (k = 0; k < np; k++)
   {
      if ((i = j = map[k]) < 0) continue;

#ifndef COEFICIENTES_CONSTANTES
      // Genera los coeficientes variables
      genCoefVar(k);
#endif

      m = ntype[k] >> 12;    // node multiplicity
      A(i, j) = coef[0] / m;
      j = (k + 1 < np ? map[k + 1] : -1);
      if (j >= 0 && ((k + 1) % Mx) != 0)
      {
         m1 = (((m1 = (ntype[k + 1] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k+1)
         A(i, j) = coef[1] / m1;
         b1 = (i > j ? i - j : j - i);
         if (b1 > b2) b2 = b1;
      }
      j = (k - 1 >= 0 ? map[k - 1] : -1);
      if (j >= 0 && (k % Mx) != 0)
      {
         m1 = (((m1 = (ntype[k - 1] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k-1)
         A(i, j) = coef[2] / m1;
         b1 = (i > j ? i - j : j - i);
         if (b1 > b2) b2 = b1;
      }
      if (nDim > 1)
      {
         j = (k + Mx < np ? map[k + Mx] : -1);
         if (j >= 0 && (nDim < 3 || (k + Mx < (k / Mxy)*Mxy + Mxy)))
         {
            m1 = (((m1 = (ntype[k + Mx] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k+Mx)
            A(i, j) = coef[3] / m1;
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
         j = (k - Mx >= 0 ? map[k - Mx] : -1);
         if (j >= 0 && (nDim < 3 || (k - Mx >= (k / Mxy)*Mxy)))
         {
            m1 = (((m1 = (ntype[k - Mx] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k-Mx)
            A(i, j) = coef[4] / m1;
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
      }
      if (nDim > 2)
      {
         j = (k + Mxy < np ? map[k + Mxy] : -1);
         if (j >= 0)
         {
            m1 = (((m1 = (ntype[k + Mxy] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k+Mxy)
            A(i, j) = coef[5] / m1;
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
         j = (k - Mxy >= 0 ? map[k - Mxy] : -1);
         if (j >= 0)
         {
            m1 = (((m1 = (ntype[k - Mxy] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k-Mxy)
            A(i, j) = coef[6] / m1;
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
      }
   }




   if (fac != 0.0) 
   {
   
      int i, j, k, eq1 = 0, eq2 = 0;
      for (i = 0; i < np; i++)
      {
         if ((eq1 = map[i]) < 0) continue;
         for (j = 0; j < np; j++)
         {
            if ((eq2 = map[j]) < 0) continue;
            if (isDual(i) && isDual(j)) A(eq1,eq2) += fac;
            k = (i > j ? i - j : j - i);
            if (k > b2) b2 = k;
         }
      }
   }
      

   if (ind == 0) {
    gmm::copy(A,MapInt);
    if (bsym) {
#ifdef PRECONDICIONADOR_ILDLTT
       gmm::ildltt_precond<gmm::csr_matrix<double> > PR(MapInt, 10, 10e-4);
       S_PR_MapInt = PR;
#endif
    } else {
#ifdef PRECONDICIONADOR_ILDLTT
      gmm::ildltt_precond<gmm::csr_matrix<double> > PR(MapInt, 10, 10e-4);
      NS_PR_MapInt = PR;
#endif
    }
//~ gmm::diagonal_precond<gmm::csr_matrix<double> > PR(MapInt);
//~ gmm::ildltt_precond<gmm::csr_matrix<double> > PR(MapInt, 10, 10e-4);
//~ gmm::ilut_precond<gmm::csr_matrix<double> > PR(MapInt, 10, 10e-4);
//~ gmm::ilutp_precond<gmm::csr_matrix<double> > PR(MapInt, 10, 10e-4);
      TamSis0 = n;
   } else {
      gmm::copy(A,MapFull);
      if (bsym) {
#ifdef PRECONDICIONADOR_ILDLTT
         gmm::ildltt_precond<gmm::csr_matrix<double> > PR(MapFull, 10, 10e-4);
         S_PR_MapFull = PR;
#endif
      } else {
#ifdef PRECONDICIONADOR_ILDLTT
         gmm::ildltt_precond<gmm::csr_matrix<double> > PR(MapFull, 10, 10e-4);
         NS_PR_MapFull = PR;
#endif
      }
      
//~ gmm::diagonal_precond<gmm::csr_matrix<double> > PR(MapFull);
//~ gmm::ildltt_precond<gmm::csr_matrix<double> > PR(MapFull, 10, 10e-4);
//~ gmm::ilut_precond<gmm::csr_matrix<double> > PR(MapFull, 10, 10e-4);
//~ gmm::ilutp_precond<gmm::csr_matrix<double> > PR(MapFull, 10, 10e-4);
      TamSis1 = n;
   }
   fflush(stdout);
}

#else
Solvable *RectSub::genInverse(int *map, ldouble fac)
{
   ce.nameClassFunct("RectSub", "genInverse");

   int b2 = 0, b1 = 0, Mx = M[0] + 1, Mxy = 0, i = 0, j = 0, k = 0, n = 0, m = 0, m1 = 0;
   Mxy = (nDim > 1 ? (M[0] + 1) * (M[1] + 1) : 0);
   for (int r = 0; r < np; r++) if (map[r] >= 0) n++;

   MatrizDispersa *A;

   int ban = nDim == 2 ? 5 : 7;
   A = new (nothrow) MatrizDispersa(n, n, ban, "A");
   if (A == NULL) ce.memoryError("A");

   for (k = 0; k < np; k++)
   {
      if ((i = j = map[k]) < 0) continue;

#ifndef COEFICIENTES_CONSTANTES
      // Genera los coeficientes variables
      genCoefVar(k);
#endif

      m = ntype[k] >> 12;    // node multiplicity
      A->asigna(i, j, coef[0] / m);
      j = (k + 1 < np ? map[k + 1] : -1);
      if (j >= 0 && ((k + 1) % Mx) != 0)
      {
         m1 = (((m1 = (ntype[k + 1] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k+1)
         A->asigna(i, j, coef[1] / m1);
         b1 = (i > j ? i - j : j - i);
         if (b1 > b2) b2 = b1;
      }
      j = (k - 1 >= 0 ? map[k - 1] : -1);
      if (j >= 0 && (k % Mx) != 0)
      {
         m1 = (((m1 = (ntype[k - 1] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k-1)
         A->asigna(i, j, coef[2] / m1);
         b1 = (i > j ? i - j : j - i);
         if (b1 > b2) b2 = b1;
      }
      if (nDim > 1)
      {
         j = (k + Mx < np ? map[k + Mx] : -1);
         if (j >= 0 && (nDim < 3 || (k + Mx < (k / Mxy)*Mxy + Mxy)))
         {
            m1 = (((m1 = (ntype[k + Mx] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k+Mx)
            A->asigna(i, j, coef[3] / m1);
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
         j = (k - Mx >= 0 ? map[k - Mx] : -1);
         if (j >= 0 && (nDim < 3 || (k - Mx >= (k / Mxy)*Mxy)))
         {
            m1 = (((m1 = (ntype[k - Mx] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k-Mx)
            A->asigna(i, j, coef[4] / m1);
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
      }
      if (nDim > 2)
      {
         j = (k + Mxy < np ? map[k + Mxy] : -1);
         if (j >= 0)
         {
            m1 = (((m1 = (ntype[k + Mxy] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k+Mxy)
            A->asigna(i, j, coef[5] / m1);
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
         j = (k - Mxy >= 0 ? map[k - Mxy] : -1);
         if (j >= 0)
         {
            m1 = (((m1 = (ntype[k - Mxy] >> 12)) > m) ? m : m1);  // m1 = min(mult)(k, k-Mxy)
            A->asigna(i, j, coef[6] / m1);
            b1 = (i > j ? i - j : j - i);
            if (b1 > b2) b2 = b1;
         }
      }
   }
   // Visualiza la matriz completa (No es como se almacena)
   //A->visualiza(1);

   Solvable *ss;

   if (fac != 0.0)  {
      int i, j, k, eq1 = 0, eq2 = 0;
      for (i = 0; i < np; i++)
      {
         if ((eq1 = map[i]) < 0) continue;
         for (j = 0; j < np; j++)
         {
            if ((eq2 = map[j]) < 0) continue;
            if (isDual(i) && isDual(j)) A->asigna(eq1, eq2, A->retorna(eq1, eq2) + fac);
            k = (i > j ? i - j : j - i);
            if (k > b2) b2 = k;
         }
      }
   }
#ifdef USAR_LOCAL_MET_ITER
   // Se usan metodos iterativos para ahorar memoria
   if (bsym)
   {
      // Resolucion del sistema lineal por medio del metodo iterativo CGM
      ss = new (nothrow) ICGM(n, A, EPSILON_LOCAL, NMAXITER_LOCAL);
      if (ss == NULL) ce.memoryError("ss in ICGM");
   }
   else
   {
      // Resolucion del sistema lineal por medio del metodo GMRES
      ss = new (nothrow) IDQGMRES(n, A, (int) sqrt(n), EPSILON_LOCAL, NMAXITER_LOCAL);
      if (ss == NULL) ce.memoryError("ss in IDQGMRES");
   }
#else
   // Metodos directos son mas rapidos en matrice pequeñas a medianas
   //~ // Resolucion del sistema lineal por medio del metodo BanSolve
   ss = new (nothrow) BandSolve(n, A);
   if (ss == NULL) ce.memoryError("ss in BandSolve");
   // Visualiza la matriz completa (No es como se almacena)
   //~ A->visualiza(1);
   
   delete A;
#endif
   return ss;

}
#endif


// Genera el numero de coordenada dentro del subdominio
void RectSub::genNcoord(int n, int *coord, int *N)
{
   for (int i = nDim - 1; i >= 0; i--)
   {
      coord[i] = (i == 0 ? n : n / N[i - 1]);
      //~ printf("\n==> %d  %d  %d",coord[i], i, N[i-1]);
      if (i > 0) n -= coord[i] * N[i - 1];
   }
}


// Genera la coordenada por subdominio
void RectSub::getCoord(int m, ldouble *x)
{
   int i;
   genNcoord(m, coord, M1);
   for (i = 0; i < nDim; i++)
   {
      x[i] = domain[i][0] + coord[i] * (domain[i][1] - domain[i][0]) / (ldouble)M[i];
      //~ printf("\n==> %d  %d %lf %d %lf %lf %d", m, i, x[i], coord[i], domain[i][0], domain[i][1], M[i]);
   }
}


vector<InternalBd*> RectSub::getInternalBd(void)
{
   ce.nameClassFunct("RectSub", "getInternalBd");

   vector<InternalBd*> intBd;
   intBd.resize(nBd);

   int nbd = 0, st;
   for (int m = 0; m < np; m++)
   {
      if ((ntype[m] & INTBD) == 0) continue;
      ldouble *y = new (nothrow) ldouble[nDim];
      if (y == NULL) ce.memoryError("y");
      getCoord(m, y);
      st = ((ntype[m] & PRIMAL) != 0 ? -1 : 0);
      intBd[nbd] = new (nothrow) InternalBd(id, m, nbd, st, nDim, y);
      if (intBd[nbd] == NULL) ce.memoryError("pintBd[nbd]", nbd);
      nbd++;
   }

   return intBd;
}


void RectSub::getPrimals(int sc, ldouble *u)
{
   for (int i = 0; i < nBd; i++)
      if ((ntype[bdMap[i]] & PRIMAL) != 0) u[i] = scr[sc][bdMap[i]];
}


void RectSub::getValues(int sc, ldouble *u)
{
   for (int i = 0; i < nBd; i++)
      if ((ntype[bdMap[i]] & DUAL) != 0) u[i] = scr[sc][bdMap[i]];
}


#ifdef GMM
void RectSub::inverse(int sp, int sc1, int sc2)
{   
   /* Computes scr[sc2][] = A(sp)-1*scr[sc1][].
    Erases other elements from scr[sc2][]  */
   int TamSis;
   if (sp == 0) TamSis = TamSis0;
    else TamSis = TamSis1;

   gmm::identity_matrix PS;   // Optional scalar product for cg
   gmm::iteration iter(EPSILON_LOCAL);// Iteration object with the max residu
   iter.set_maxiter(NMAXITER_LOCAL);
   size_t restart = (int) sqrt(TamSis);       // restart parameter for GMRES
   

   // Vectores
   std::vector<double> x(TamSis), b(TamSis);
   
   int i, j;
   int *map = (sp == 0 ? mapInt : mapFull);
   for (i = 0; i < np; i++)
      if ((j = map[i]) >= 0) b[j] = scr[sc1][i];
#ifdef SIN_PRECONDICIONADOR
   if (sp == 0) {
      if (bsym) gmm::cg(MapInt, x, b, PS, PR, iter);
       else gmm::gmres(MapInt, x, b, PR, restart, iter);
   } else {
      if (bsym) gmm::cg(MapFull, x, b, PS, PR, iter);
       else gmm::gmres(MapFull, x, b, PR, restart, iter);
   }
#endif
#ifdef PRECONDICIONADOR_ILDLTT
   if (sp == 0) {
      if (bsym) gmm::cg(MapInt, x, b, PS, S_PR_MapInt, iter);
       else gmm::gmres(MapInt, x, b, NS_PR_MapInt, restart, iter);
   } else {
      if (bsym) gmm::cg(MapFull, x, b, PS, S_PR_MapFull, iter);
       else gmm::gmres(MapFull, x, b, NS_PR_MapFull, restart, iter);
   }
#endif   
   
//printf("\nIteraciones Locales %d",iter.get_iteration());
   for (i = 0; i < np; i++)
      scr[sc2][i] = ((j = map[i]) >= 0 ? x[j] : 0.0);
}


#else
void RectSub::inverse(int sp, int sc1, int sc2)
{
   /* Computes scr[sc2][] = A(sp)-1*scr[sc1][].
    Erases other elements from scr[sc2][]  */
   int i, j;
   int *map = (sp == 0 ? mapInt : mapFull);
   for (i = 0; i < np; i++)
      if ((j = map[i]) >= 0) Y[j] = scr[sc1][i];
   inv[sp]->solve(X, Y);
   for (i = 0; i < np; i++)
      scr[sc2][i] = ((j = map[i]) >= 0 ? X[j] : 0.0);
}
#endif

// Genera el tipo de nodo de cada nodo del subdominio
void RectSub::genNtype(Primal &primal)
{
   //printf("\n\ngenNtype CoordN % d  %d  %d, M1 %d %d",coordN[0], coordN[1], coordN[2], M1[0], M1[1]);
   int m, type = 0, r = 0;
   for (m = 0; m < np; m++)
   {
      genNcoord(m, coord, M1);
      //printf("\nind = %d, Coord % d  %d %d",m,coord[0], coord[1], coord[2]);
      type = 0;
      if (isKnown(coord)) type |= KNOWN;
      else if (isInterior(coord)) type |= INTERIOR;
      else
      {
         isIntBd(coord);
         type |= INTBD;
         r = nodeType(coord);
         if (r == nDim) type |= VERTEX;
         else if (nDim == 3 && r == 1) type |= FACE;
         else type |= EDGE;
      }
      if ((type & INTBD) != 0 && primal.isPrimal(type, coordN, coord))
         type |= PRIMAL;
      if ((type & INTBD) != 0 && (type & PRIMAL) == 0) type |= DUAL;
      if ((type & INTERIOR) != 0) type |= (1 << 12);  // multiplicity of interior nodes
      ntype[m] = type;
      //printf("\nNtype  %d",type);
   }
}


// Es nodo conocido
bool RectSub::isKnown(int *coord)
{
   int i, sw = 0;
   if (nDim != 3)
   {
      int ndIniSub, nd, ndx, ndy , ndz, ndxy, ndxyz;
      ndx = mesh[0] * mesh[3] + 1;
      ndy = mesh[1] * mesh[4] + 1;
      //~ ndz = mesh[2]*mesh[5]+1;
      ndxy = ndx * ndy;
      ndxyz = ndxy * mesh[2] * mesh[5] + 1;

      ndIniSub = (ndxy * coordN[2]) + (ndx * mesh[4] * coordN[1]) + (mesh[3] * coordN[0]);
      nd = ndIniSub + (ndxyz * coord[2]) +  (ndx * coord[1]) + coord[0];
      //printf("\n CoordN  %d %d %d  Base %d Nodo %d ",coordN[0], coordN[1], coordN[2], ndIniSub, nd);



      if(nd <= ndx) sw = 1;
      if(nd >= (ndx * (mesh[1]*mesh[4])) ) sw = 1;
      if((nd % ndx) == 0) sw = 1;
      if(((nd + 1) % ndx) == 0) sw = 1;

      if(sw == 1)
      {
         //printf ("\nConocido nod = %d (%d, %d)\n",nd,coord[0], coord[1]);
         return true;
      }
   }
   else
   {
      for (i = 0; i < nDim; i++)
      {
         if ((coord[i] == 0 && coordN[i] == 0) || (coord[i] == M[i] && coordN[i] == M[i] - 1))
         {
            return true;
         }
      }
   }
   return false;
}

// Es nodo interior
bool RectSub::isInterior(int *coord)
{
   int i, sw = 0, ndl;
   if (nDim != 3)
   {
      int bas, ndx, nd;
      ndx = mesh[0] * mesh[3] + 1;
      bas = (mesh[3] * coordN[0]) + (ndx * mesh[4]) * coordN[1];
      nd = bas + ndx * coord[1] + coord[0];
      //printf("\n CoordN  %d %d Base %d Nodo %d ",coordN[0], coordN[1], bas, nd);

      if(nd <= ndx) sw = 1;
      if(nd >= (ndx * (mesh[1]*mesh[4])) ) sw = 1;
      if((nd % ndx) == 0) sw = 1;
      if(((nd + 1) % ndx) == 0) sw = 1;


      if(sw != 1)
      {
         ndl = (mesh[3] + 1) * coord[1] + coord[0];

         if (ndl <= (mesh[3] + 1)) return false;
         if (ndl >= ((mesh[3] + 1)*mesh[4])) return false;
         if (((ndl + 1)  % (mesh[3] + 1)) == 0 ||  (ndl  % (mesh[3] + 1)) == 0) return false;
         //printf ("\nInterior nod glob = %d, nodo Local %d  (%d, %d)\n",nd, ndl, coord[0], coord[1]);
         return true;
      }
   }
   else
   {
      for (i = 0; i < nDim; i++)
         if (coord[i] == 0 || coord[i] == M[i]) return false;
      return true;
   }
   return false;
}

// Es nodo de la frontera interior
bool RectSub::isIntBd(int *coord)
{
   int i, sw = 0, ndl;
   if (nDim != 3)
   {
      int bas, ndx, nd;
      ndx = mesh[0] * mesh[3] + 1;
      bas = (mesh[3] * coordN[0]) + (ndx * mesh[4]) * coordN[1];
      nd = bas + ndx * coord[1] + coord[0];
      //printf("\n CoordN  %d %d Base %d Nodo %d ",coordN[0], coordN[1], bas, nd);

      if(nd <= ndx) sw = 1;
      if(nd >= (ndx * (mesh[1]*mesh[4])) ) sw = 1;
      if((nd % ndx) == 0) sw = 1;
      if(((nd + 1) % ndx) == 0) sw = 1;

      if(sw != 1)
      {
         ndl = (mesh[3] + 1) * coord[1] + coord[0];
         sw = 0;
         if (ndl <= (mesh[3] + 1)) sw = 1;
         if (ndl >= ((mesh[3] + 1)*mesh[4])) sw = 1;
         if (((ndl + 1)  % (mesh[3] + 1)) == 0 ||  (ndl  % (mesh[3] + 1)) == 0) sw = 1;

         if (sw == 1)
         {
            //printf ("\nFrontera Interior nod glob = %d, nodo Local %d  (%d, %d)\n",nd, ndl, coord[0], coord[1]);
            return true;
         }
      }
   }
   else
   {
      for (i = 0; i < nDim; i++)
         if (coord[i] == 0 || coord[i] == M[i])
            if (coordN[i] != 0 && coordN[i] != N[i]) return true;
   }
   return false;
}

// Revisa el tipo de nodo (Vertice, Cara, Arista)
int RectSub::nodeType(int *coord)
{

   int i, sw = 0, ndl;
   if(nDim != 3)
   {
      int bas, ndx, nd;
      ndx = mesh[0] * mesh[3] + 1;
      bas = (mesh[3] * coordN[0]) + (ndx * mesh[4]) * coordN[1];
      nd = bas + ndx * coord[1] + coord[0];
      //printf("\n CoordN  %d %d Base %d Nodo %d ",coordN[0], coordN[1], bas, nd);

      if(nd <= ndx) sw = 1;
      if(nd >= (ndx * (mesh[1]*mesh[4])) ) sw = 1;
      if((nd % ndx) == 0) sw = 1;
      if(((nd + 1) % ndx) == 0) sw = 1;

      if(sw != 1)
      {
         ndl = (mesh[3] + 1) * coord[1] + coord[0];
         sw = 0;
         if (ndl <= (mesh[3] + 1)) sw = 1;
         if (ndl >= (mesh[3] * (mesh[4] + 1))) sw = 1;
         if (((ndl + 1)  % (mesh[3] + 1)) == 0 ||  (ndl  % (mesh[3] + 1)) == 0) sw = 1;

         if (sw == 1)
         {
            if (ndl == 0 || ndl == mesh[3] || ndl == (mesh[3] + 1)*mesh[4] || ndl == ((mesh[3] + 1)*mesh[4] + mesh[3]))
            {
               //printf ("\nVertice nod glob = %d, nodo Local %d  (%d, %d)\n",nd, ndl, coord[0], coord[1]);
               return 2;
            }
         }
      }
   }
   else
   {
      int n = 0;
      for (int i = 0; i < nDim; i++)
      {
         if (coord[i] == 0 || coord[i] == M[i]) n++;
         //printf("\n(%d)  Coord %d  M %d  n %d",i, coord[i], M[i], n);
      }
      return n;
   }
   return 0;
}

void RectSub::knownValues(int s1)
{
   FunctionV *g = op->getG();
   for (int i = 0; i < np; i++)
   {
      if (!isKnown(i)) continue;
      getCoord(i, x);
      scr[s1][i] = g->eval(0, x);
   }
}

void RectSub::multOp(int s1, int s2)
{
   int i, Mx = M[0] + 1, Mxy, m = 0, m1 = 0;
   Mxy = (nDim > 2 ? (M[0] + 1) * (M[1] + 1) : 0);
   for (i = 0; i < np; i++)
   {
      if (isKnown(i))
      {
         scr[3][i] = 0.0;
         continue;
      }

#ifndef COEFICIENTES_CONSTANTES
      // Genera los coeficientes variables
      genCoefVar(i);
#endif

      m = ntype[i] >> 12;
      scr[3][i] = coef[0] * scr[s1][i] / m;
      if (i + 1 < np  && ((i + 1) % Mx) != 0)
      {
         m1 = ((m1 = (ntype[i + 1] >> 12)) > m ? m : (m1 == 0 ? 1 : m1));
         scr[3][i] += coef[1] * scr[s1][i + 1] / m1;
      }
      if (i > 0 && (i % Mx) != 0)
      {
         m1 = ((m1 = (ntype[i - 1] >> 12)) > m ? m : (m1 == 0 ? 1 : m1));
         scr[3][i] += coef[2] * scr[s1][i - 1] / m1;
      }
      if (nDim > 1)
      {
         if (i + Mx < np && (nDim < 3 || ((i + Mx) < (i / Mxy)*Mxy + Mxy)))
         {

            m1 = ((m1 = (ntype[i + Mx] >> 12)) > m ? m : (m1 == 0 ? 1 : m1));
            scr[3][i] += coef[3] * scr[s1][i + Mx] / m1;
         }
         if (i - Mx >= 0 && (nDim < 3 || ((i - Mx) >= (i / Mxy) *Mxy)))
         {
            m1 = ((m1 = (ntype[i - Mx] >> 12)) > m ? m : (m1 == 0 ? 1 : m1));
            scr[3][i] += coef[4] * scr[s1][i - Mx] / m1;
         }
      }
      if (nDim > 2)
      {
         if (i + Mxy < np)
         {
            m1 = ((m1 = (ntype[i + Mxy] >> 12)) > m ? m : (m1 == 0 ? 1 : m1));
            scr[3][i] += coef[5] * scr[s1][i + Mxy] / m1;
         }
         if (i - Mxy >= 0)
         {
            m1 = ((m1 = (ntype[i - Mxy] >> 12)) > m ? m : (m1 == 0 ? 1 : m1));
            scr[3][i] += coef[6] * scr[s1][i - Mxy] / m1;
         }
      }
   }
   for (i = 0; i < np; i++) scr[s2][i] = scr[3][i];
}


void RectSub::printMat(const char *s, ldouble **A, int tm)
{
   int i, j;
   printf("\n%s\n", s);
   for (i = 0; i < tm; i++)
   {
      printf("%d ", i);
#ifdef __Double__
      for (j = 0; j < tm; j++) printf("%f ", A[i][j]);
#else
      for (j = 0; j < tm; j++) printf("%Lf ", A[i][j]);
#endif
      printf("\n");
   }
}

void RectSub::printMult(void)
{
   ce.nameClassFunct("RectSub", "printMult");

   int i, j;
   ldouble **A = new (nothrow) ldouble *[np];
   if (A == NULL) ce.memoryError("A");
   for (i = 0; i < np; i++)
   {
      A[i] = new (nothrow) ldouble[np];
      if (A[i] == NULL) ce.memoryError("A[i]", i);
   }
   for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) A[i][j] = 0.0;

   for (i = 0; i < np; i++)
   {
      for (j = 0; j < np; j++) scr[0][j
                                        ] = (i == j ? 1.0 : 0.0);
      multOp(0, 1);
      for (j = 0; j < np; j++) A[i][j] = scr[1][j];
   }
   printMat("MULT", A, np);

   for (i = 0; i < np; i++) delete []A[i];
   delete []A;
}



void RectSub::rhs(int sc)
{
   int i, m;
   ldouble val;
   FunctionV *f = op->getF();
   FunctionV *g = op->getG();
   f->dimension(nDim);
   g->dimension(nDim);
   for (i = 0; i < np; i++)
      if ((ntype[i] & KNOWN) != 0)
      {
         getCoord(i, x);
         scr[2][i] = (g == NULL ? 0.0 : g->eval(0, x));
      }
      else scr[2][i] = 0.0;

   multOp(2, 3);
   for (i = 0; i < np; i++) scr[sc][i] = -scr[3][i];
   for (i = 0; i < np; i++)
   {
      if ((ntype[i] & KNOWN) != 0) continue;
      getCoord(i, x);
      m = (((ntype[i] & INTBD) != 0) ? ntype[i] >> 12 : 1);
      val = (f == NULL ? scr[sc][i] / m : (hfac * f->eval(0, x) + scr[sc][i]) / m);
      scr[sc][i] = val;
   }
}

void RectSub::setPrimals(int sc, ldouble *u)
{
   for (int i = 0; i < nBd; i++)
      if ((ntype[bdMap[i]] & PRIMAL) != 0) scr[sc][bdMap[i]] = u[i];
}


void RectSub::setValues(int sc, ldouble *u)
{
   for (int i = 0; i < nBd; i++)
      if ((ntype[bdMap[i]] & DUAL) != 0) scr[sc][bdMap[i]] = u[i];
}


void RectSub::print(const char *s, int sc)
{
   // Salida por subdominio
   printf("\n%s %d", s, id);
#ifdef __Double__
   for (int i = 0; i < np; i++) printf("\n%03d %1.16e", i, scr[sc][i]);
#else
   for (int i = 0; i < np; i++) printf("\n%03d %1.16Le", i, scr[sc][i]);
#endif
   printf("\n");

   fflush(stdout);
}


void RectSub::print(int sc)
{
   // Salida en la que se muestra las coordenadas y la solución
   int i;
   ldouble *x = new (nothrow) ldouble[nDim];
   switch (nDim)
   {
   case 1:
      for (i = 0; i < np; i++)
      {
         getCoord(i, x);
#ifdef __Double__
         printf("\n %+1.16e  %+1.16e", x[0], scr[sc][i]);
#else
         printf("\n %+1.16Le  %+1.16Le", x[0], scr[sc][i]);
#endif
      }
      break;
   case 2:
      for (i = 0; i < np; i++)
      {
         getCoord(i, x);
#ifdef __Double__
         printf("\n %+1.16e %+1.16e  %+1.16e", x[0], x[1], scr[sc][i]);
#else
         printf("\n %+1.16Le %+1.16Le  %+1.16Le", x[0], x[1], scr[sc][i]);
#endif
      }
      break;
   default:
      for (i = 0; i < np; i++)
      {
         getCoord(i, x);
#ifdef __Double__
         printf("\n %+1.16e %+1.16e %+1.16e  %+1.16e", x[0], x[1], x[2], scr[sc][i]);
#else
         printf("\n %+1.16Le %+1.16Le %+1.16Le  %+1.16Le", x[0], x[1], x[2], scr[sc][i]);
#endif
      }
   };
   delete []x;

   fflush(stdout);
}
