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



#include <string>
#include <time.h>
#include "DPMethod.hpp"




void DPMethod::initialize(void)
{
   ce.nameClassFunct("DPMethod", "inicializa");

   time(&time0);


   //inter = new Interchange(*props);
   iniInterchage();

   nOmega = inter->getnOmega();
   nDim = inter->getnDim();

   genInverse(0);
   time(&time1);
   printv = props->getInt("p", 0);
   dualp = new (nothrow) DualPrimal(*inter);
   if (dualp == NULL) ce.memoryError("dualp");
   time(&time2);
   nDual = dualp->getNDual();
   rhss = new (nothrow) ldouble[nDual];
   if (rhss == NULL) ce.memoryError("rhss");
   u = new (nothrow) ldouble[nDual];
   if (u == NULL) ce.memoryError("u");
   int i;
   for (i = 0; i < nDual; i++)
   {
      rhss[i] = 0.0;
      u[i] = 0.0;
   }
}


/* Computes inverse in each subdomain or resuses an inverse for floating domains
     If type = 0, normal case, matrix = S, type = 1, then matrix = M
*/
void DPMethod::genInverse(int type)
{
   for (int e = 0; e < inter->getnOmega(); e++) inter->genInv(e, type);
}

// Imprime la solucion al problema
// (1) Visualiza la solución de nodos duales
// (2) Visualiza la solución de todos los nodos por subdominio
// (3) Visualiza las coordenadas y su solución
void DPMethod::print(ldouble *u)
{
   int e;
   ldouble *x = new (nothrow) ldouble[nDim];
   if (x == NULL) ce.memoryError("x");
   for (e = 0; e < nDim; e++) x[e] = 0.0;

   time(&time3);
   printf("\nprint dimension (%d) impresion type (%d)", nDim, printv);
   if (printv == 1)
   {
      // Visualiza la solución de nodos duales
      int i = 0;
      size_t j;
      for (j = 0; j < inter->bds->bdDuals.size(); j++)
      {
         inter->getCoordNode(inter->bds->bdDuals[j][0]->subd, inter->bds->bdDuals[j][0]->node, x);
#ifdef __Double__
         printf("\n%4d %s %+1.16e", i, prCoord(x), u[inter->bds->bdDuals[j][0]->index]);
#else
         printf("\n%4d %s %+1.16Le", i, prCoord(x), u[inter->bds->bdDuals[j][0]->index]);
#endif
         i++;
      }
   }
   if (printv == 2)
   {
      // Visualiza la solución de todos los nodos por subdominio
      inter->print("Subdomain", 0);
   }
   if (printv == 3)
   {
      // Visualiza las coordenadas y su solución
      inter->print(0);
   }


   printf("\n%s iterations =  %d\n", solver->getName(), solver->getIter());
   printTime();
   delete []x;
}


void DPMethod::printTime(void)
{
   printf("\nLocal Inverses   %f", difftime(time1, time0));
   printf("\nPrimal Matrices  %f", difftime(time2, time1));
   printf("\nIteration Time   %f", difftime(time3, time2));
   printf("\nTotal Time       %f\n\n", difftime(time3, time0));
   fflush(stdout);
}


const char *DPMethod::prCoord(ldouble *x)
{
   char txt[1000];
   const char *r;

   string s = "(";
   for (int i = 0; i < nDim; i++)
   {
#ifdef __Double__
      sprintf(txt, "%+1.16e", x[i]);
#else
      sprintf(txt, "%+1.16Le", x[i]);
#endif
      s += (string) txt + (i == nDim - 1 ? ")" : ", ");
   }
   r = s.c_str();

   return r;
}


double DPMethod::analyticSolution(double *x)
{
   return exp(x[0] + x[1]);
}

/// Calcula el numero de condicionamiento
// Autor Alberto Rosas Medina
void DPMethod::conditionalNumber(bool symetric)
{
   /*
      int i;
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //NUMERO DE CONDICIONAMIENTO Y SENO DEL ANGULO
      //Calculo del error para u-unumerica y u
      double MaxdeError=0.0;
      double MaxdeU=0.0;
      double r = 1.0;
      // int rsubd, rindex,rnode;
      double x[3];
      double SenodelAngulo=0.0;// Numero de con matrices no simetricas
      double SenodelAngulo1=0.0;// Numero de con matrices no simetricas
      double NumerodeCondNoSym=0.0;
      double NumerodeCondNoSym1=0.0;


      double *T = new (nothrow) double[nDual];
      if (T == NULL) ce.memoryError("T");
      for (i = 0; i < nDual; i++) T[i] = 0.0;

      double *T1 = new (nothrow) double[nDual];
      if (T1 == NULL) ce.memoryError("T1");
      for (i = 0; i < nDual; i++) T1[i] = 0.0;

      double *Su = new (nothrow) double[nDual];
      if (Su == NULL) ce.memoryError("Su");
      for (i = 0; i < nDual; i++) Su[i] = 0.0;

      double *error1 = new (nothrow) double[nDual];
      if (error1 == NULL) ce.memoryError("error");
      for (i = 0; i < nDual; i++) error1[i] = 0.0;

      double *Resta = new (nothrow) double[nDual];
      if (Resta == NULL) ce.memoryError("T");
      for (i = 0; i < nDual; i++) Resta[i] = 0.0;

      for ( i = 0; i < inter->bdDuals.size(); i++)
         //for (int i = 0; i < nDual; i++)
      {
         inter->getCoordNode(inter->bdDuals[i][0]->rsubd(),inter->bdDuals[i][0]->rnode(), x);

         if(nDim==2)
         {
            //Su[inter->bdDuals[i][0]->rindex()]=sin(M_PI*r*x[0])*sin(M_PI*r*x[1]);
            //Su[inter->bdDuals[i][0]->rindex()]=sin(M_PI*x[0])*sin(M_PI*x[1]);
            //Su[inter->bdDuals[i][0]->rindex()]=exp(x[0]*x[1]);
            //Su[inter->bdDuals[i][0]->rindex()]=exp(3.0*x[0]+3.0*x[1]);
            Su[inter->bdDuals[i][0]->rindex()]=analyticSolution(x);
            // Su[inter->bdDuals[i][0]->rindex()]=sin(M_PI*x[0])*cos(M_PI*x[1]);

            //  printf("\n Solución exacta ...*......*.......*........*.........*.......*.......*...... %e , %e:\n",  x[0],x[1]);
            //printf("\n Solución exacta ...*......*.......*........*.........*.......*.......*...... %e :\n",  Su[inter->inter->bdDuals[i][0]->rindex()]);
            error1[inter->bdDuals[i][0]->rindex()]=Su[inter->bdDuals[i][0]->rindex()]-u[inter->bdDuals[i][0]->rindex()];
            // printf("\n errores ..................................................... %e :\n",  error1[inter->bdDuals[i][0]->rindex()]);
            // printf("\n print(Su) .....................................................  :\n");
            // print(Su);
         }
         else
         {
            //Su[bdDuals[i][0]->rindex()]=sin(M_PI*r*x[0])*sin(M_PI*r*x[1])*sin(M_PI*r*x[2]);
            Su[inter->bdDuals[i][0]->rindex()]=exp(x[0]+x[1]+x[2]);
            // Su[bdDuals[i][0]->rindex()]=exp(x[0]*x[1]*x[2]);
            error1[inter->bdDuals[i][0]->rindex()]=Su[inter->bdDuals[i][0]->rindex()]-u[inter->bdDuals[i][0]->rindex()];
            //NO error[i]=(u[bdDuals[i][0]->rindex()]-Su[i]);
         }
         if (fabs(error1[inter->bdDuals[i][0]->rindex()])>MaxdeError) MaxdeError=fabs(error1[inter->bdDuals[i][0]->rindex()]);


         if (fabs(Su[inter->bdDuals[i][0]->rindex()])>MaxdeU) MaxdeU=fabs(Su[inter->bdDuals[i][0]->rindex()]);
      }
   // printf("\n Error u numerica - u exacta_______________________________________NUEVO \n");
   //print(error1);
   //~ printf("\n Solucion analitica \n");
   //~ printf("\n NDualSize %d :\n", inter->bdDuals.size());
   //~ printf("\n NDual %d :\n", nDual);
   //~ print(Su);
   //printf("\n Error u numerica - u exacta \n");
   //print(error);
   //printf("\n Maximo del Error es %e\n",MaxdeError);




   // /////////////////////////////////////////////////////////////////////
   //Norma A  del Vector error=u-uk
      double val = 0.0;
      double val1 = 0.0;
      double val2 = 0.0;
      double val3 = 0.0;
      double Valerror=0.0;
      double ValSu=0.0;
      double valornormacuadrada=0.0;
      double valeRROR=0.0;
      double NormaEuclidianaDelError =0.0;
      double NormaEuclidianaDeU =0.0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //NUEVOS CALCULOS DE LAS NORMAS
      //multOp(error1, T);
      ///////////////////////////////////////////////////////////////////////////
      // CALCULO DE  Uk A Uk
      multOp(u, T);
   // printf("\n vector de A x error \n");
   //print(T);
      // CALCULO DE  Uk A Uk
      for (int i = 0; i < inter->bdDuals.size(); i++) val += T[inter->bdDuals[i][0]->rindex()]*u[inter->bdDuals[i][0]->rindex()];
   /////////////////////////////////////
   // CALCULO DE  U A Uk
      for (int i = 0; i < inter->bdDuals.size(); i++) val1 += T[inter->bdDuals[i][0]->rindex()]*Su[inter->bdDuals[i][0]->rindex()];

   ///////////////////////////////////////////////////////////////////////////
   // CALCULO DE   Su A Su
      multOp(Su, T1);
      for (int i = 0; i < inter->bdDuals.size(); i++) val2 += T1[inter->bdDuals[i][0]->rindex()]*Su[inter->bdDuals[i][0]->rindex()];
      // CALCULO DE Uk A Su
      for (int i = 0; i < inter->bdDuals.size(); i++) val3 += T1[inter->bdDuals[i][0]->rindex()]*u[inter->bdDuals[i][0]->rindex()];

   // CALCULO DE ||U||
      ValSu = sqrt(val2);
      printf("\n Maximo del Error es %e\n",MaxdeError);// ERROR MÁXIMO EN NORMA INFINITO
      printf("\n Norma A  cuadrada de la Solución analitica es %1.16e: \n", val2);
      printf("\n Norma A de la solucion analitica es %1.16e: \n", ValSu);

   //          CALCULO DE   UAU-UkAU-UAUk+UkAUk
      printf("\n valores de uAu\n");
      printf("\n val2 %1.16e: \n", val2);
      printf("\n val3 %1.16e: \n", val3);
      printf("\n val1 %1.16e: \n", val1);
      printf("\n val %1.16e: \n", val);
   //valornormacuadrada= val2-val3-val1+val;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //Calculo de U-Uk A U-Uk
      multOp(error1, Resta);
      for (int i = 0; i < inter->bdDuals.size(); i++) valeRROR += error1[inter->bdDuals[i][0]->rindex()]*Resta[inter->bdDuals[i][0]->rindex()];
      printf("\n Valerror%1.16e: \n", valeRROR);
   //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   // calculo drecto de U-Uk A U-Uk
      valornormacuadrada= valeRROR;

   // CALCULO DE LA NORMA DE || U-Uk||
      Valerror=sqrt(valornormacuadrada);//
      printf("\n Norma A cuadrada del error es %1.16e: \n", valornormacuadrada);
      printf("\n Norma A del error es %1.16e: \n", Valerror);
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // Calculo del valor de q
      int NIter=0;
      NIter=solver->getIter();
      double NumerodeCond=0.0;
      double q=0.0;
      double tq=0.0;
      double tITer=0.0;
      double TemSeno=0.0;
      TemSeno=1.0/(double)NIter;

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // CAULCULO DE LA NORMA DEL ERROR CON LA NORMA 2
      for (int i = 0; i < inter->bdDuals.size(); i++)  NormaEuclidianaDelError += error1[inter->bdDuals[i][0]->rindex()]*error1[inter->bdDuals[i][0]->rindex()];

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // CAULCULO DE LA NORMA DE LA SOLUCION ANALITICA CON LA NORMA 2
      for (int i = 0; i < inter->bdDuals.size(); i++)  NormaEuclidianaDeU += Su[inter->bdDuals[i][0]->rindex()]*Su[inter->bdDuals[i][0]->rindex()];

      double RaizdelaNormaDelError=0.0;
      double RaizdelaNormaDeU=0.0;

      RaizdelaNormaDelError=pow(NormaEuclidianaDelError,0.5);
      RaizdelaNormaDeU=pow(NormaEuclidianaDeU,0.5);

      if (symetric)
      {
         //tq=(Valerror/ValSu); BUENO
         tq=(Valerror/(2*ValSu));
         printf("\n El numero de tq es %1.16e: \n", tq);
         //tITer=1.0/(double)(2.0*NIter); BUENO
         tITer=1.0/NIter;

         q=pow(tq,tITer);
         NumerodeCond=pow( ((1.0+q)/(1.0-q)),2.0);
         printf("\n El numero de q es %1.16e: \n", q);
         printf("\n El numero de condicionamiento es %1.16e: \n", NumerodeCond);
         printf("\n El numero de ITER es %d: \n", NIter);
         printf("\n La raiz  es %e: \n", tITer);
         //printf("\n la division es es %1.16e: \n", tq);


      }

      else
      {

         SenodelAngulo=pow((MaxdeError/MaxdeU),TemSeno);///////NORMA INFINITO
         SenodelAngulo1=pow((RaizdelaNormaDelError/RaizdelaNormaDeU),TemSeno);
         NumerodeCondNoSym1=pow((SenodelAngulo1+1)/(1-SenodelAngulo1),2.0);
         NumerodeCondNoSym=pow((SenodelAngulo+1)/(1-SenodelAngulo),2.0);

   //printf("\n Norma infinita del error es %e \n", MaxdeError);
   // printf("\n Norma infinita del vector solucion es %e \n", MaxdeU);
         printf("\n el exponente para el seno del angulo alpha es %e \n", TemSeno);
         printf("\n seno del angulo alpha es %e \n", SenodelAngulo);
         printf("\n seno del angulo alpha efectivo es %e \n", SenodelAngulo1);
         printf("\n El numero de Condicionamiento  para no simetricas es----- %e \n", NumerodeCondNoSym1);
         printf("\n El numero de Condicionamiento es -----prmimero- %e \n", NumerodeCondNoSym);
      }
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   */
}
