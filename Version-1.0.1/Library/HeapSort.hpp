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



#ifndef __HeapSort__
#define __HeapSort__

#include <vector>
#include "InternalBd.hpp"


class HeapSort
{
   /* A heap is an array r[0:n]; however, we use only the elements r[1:n]. The elements
      r[i] satisfy the relationship r[i] >= r[2*i], r[2*i + 1]
   */
private:

   vector<InternalBd*>  r;
   int n;

public:

   HeapSort(vector<InternalBd*>  &a, int n)
   {
      r = a;
      this->n = n;
   }

   ~HeapSort()
   {
      size_t i;
      for (i = 0; i < r.size(); i++)
      {
         InternalBd *abc = r[i];
         //delete abc;
         abc = NULL;
      }
      r.clear();
   }


   inline void genHeap(void)
   {
      for (int i = (n >> 1); i >= 1; i--) siftup(i, n);
   }

   void siftup(int i, int n)
   {
      /* Creates a subheap at i assuming there are already
         subheaps at 2*i and 2*i + 1
      */

      int j;
      while ((j = (i << 1)) <= n)
      {
         if (j < n) if (r[j]->compareTo(r[j + 1]) < 0) j++;
         if (r[i]->compareTo(r[j]) < 0)
         {
            swap(i, j);
            i = j;
         }
         else break;
      }
   }

   void sort(void)
   {
      /* sorts the array r[] 1:n by first creating a heap and then successively
         taking the largest element and putting it at the current "end" (which goes from n to 2)
      */
      int i;
      genHeap();
      swap(1, n);
      for (i = n - 1; i >= 2; i--)
      {
         siftup(1, i);
         swap(1, i);
      }

   }

   inline void swap(int i, int j)
   {

      InternalBd *t = r[i];
      r[i] = r[j];
      r[j] = t;
      //~ printf("\n swap %d %d",i,j);
   }


   // Retorna el InternalBd de dato del que se trate el HeapSort
   inline InternalBd* rr(int i)
   {
      return r[i];
   }

};

#endif





