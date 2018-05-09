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



#ifndef __CreateBdNodes__
#define __CreateBdNodes__


#include "RectSub.hpp"
#include "BdNode.hpp"
#include "HeapSort.hpp"
#include "Interchange.hpp"
#include "ErrorControl.hpp"


/* This class creates the arrays BdNode[][] bdDuals, bdPrimals required by the
   class DualPrimal. A single object of this class should be created.
   The bdDuals, bdPrimals arrays should be extracted and then the object ahould
   be destroyed to return space to the system.
*/

class CreateBdNodes
{


protected:

   /// Control de errores
   ErrorControl ce;


public:

   int nD;   // number of duals
   int nP;   // number of primals
   int nDual;
   int nPrimal;
   int maxBd;    // total number of internal boundary elements from all subdomains
   int ibd;
   int ibdAll;  // total number of doundary element sets
   int *dualMult;

   vector<vector<BdNode*> > bdAll;
   vector<vector<BdNode*> > bdDuals;
   vector<vector<BdNode*> > bdPrimals;
   vector<InternalBd*>      hbd;


   CreateBdNodes(void)
   {
      dualMult = NULL;
      ibd = 0;
      ibdAll = 0;
      nD = nP = 0;
   }


   ~CreateBdNodes()
   {
      size_t k, j, i;
      for (k = 0; k < bdDuals.size(); k++)
      {
         for (j = 0; j < bdDuals[k].size(); j++)
         {
            BdNode *abc = bdDuals[k][j];
            //~ delete abc;
            abc = NULL;
         }
         bdDuals[k].clear();
      }
      bdDuals.clear();

      for (k = 0; k < bdPrimals.size(); k++)
      {
         for (j = 0; j < bdPrimals[k].size(); j++)
         {
            BdNode *abc = bdPrimals[k][j];
            //~ delete abc;
            abc = NULL;
         }
         bdPrimals[k].clear();
      }
      bdPrimals.clear();


      for (k = 0; k < bdAll.size(); k++)
      {
         for (j = 0; j < bdAll[k].size(); j++)
         {
            BdNode *abc = bdAll[k][j];
            delete abc;
            abc = NULL;
         }
         bdAll[k].clear();
      }
      bdAll.clear();


      for (i = 0; i < hbd.size(); i++)
      {
         InternalBd *abc = hbd[i];
         delete abc;
         abc = NULL;
      }
      hbd.clear();


      delete []dualMult;
      dualMult = NULL;
   }



   //~ inline int *getDualMult()
   //~ {
   //~ return dualMult;
   //~ }

   //~ inline int getMaxBd()
   //~ {
   //~ return maxBd;
   //~ }

   //~ inline int getND()
   //~ {
   //~ return nD;
   //~ }

   //~ inline int getNDuals()
   //~ {
   //~ return nDual;
   //~ }

   //~ inline int getNP()
   //~ {
   //~ return nP;
   //~ }

   //~ inline int getNPrimals()
   //~ {
   //~ return nPrimal;
   //~ }


};

#endif
