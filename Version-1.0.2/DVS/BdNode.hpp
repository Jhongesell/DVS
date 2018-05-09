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


#ifndef __BdNode__
#define __BdNode__



class BdNode
{

public:

   int subd;     // # subdomain
   int node;     //  node # in the subdomain
   int index;    // index number in DualPrimal
   int mult;     // multiplicity

   BdNode(int s, int n, int i, int m)
   {
      subd = s;
      node = n;
      index = i;
      mult = m;
   }

};

#endif
