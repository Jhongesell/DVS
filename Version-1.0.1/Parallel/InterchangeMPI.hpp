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



#ifndef  __InterchangeMPI__
#define  __InterchangeMPI__

#include "Definiciones.hpp"
#include "Interchange.hpp"
#include "EsquemaMEMPI.hpp"



/* Controls all communication and data exchange between master and subdomains */
class InterchangeMPI : public Interchange
{

private:

   /// Numero de esclavo en el que estara la tarea
   int xnp;
   /// Numero de tarea dentro del esclavo
   int indl;
   /// Arreglo para recibir mensajes
   int msa[10];
   /// Arreglo para enviar mensajes
   int mss[10];

   /// Puntero al esquema Maestro-Esclavo
   EsquemaMEMPI *ME;


public:

   /// Constructor
   InterchangeMPI(PropDef &props, EsquemaMEMPI &me);

   int getMaxBdSize(void);

   int *getNtype(int e);

   // Se agrega para poder actualizar Ntype en paralelo (ya que no se pasan datos por referencia)
   void setNtype(int e, int * arr);

   vector<InternalBd*> getInternalBd(int e);


   void calcula(int e, int node, int sp);


   /// Clear scr[sc][] in all subdomains
   void clear(int sc);

   ldouble getValue(int e, int scr1, int scr2, int node);

   /// scr[sc3][] = scr[sc1][] - scr[sc2][] in all subdomains
   void diff(int sc3, int sc1, int sc2);

   /// scr[sc2][] = A(sp)-1(scr[sc1][])
   void inverse(int sp, int sc1, int sc2);

   /// scr[sc][] = Dirichlet boundary values of all subdomains
   void knownValues(int sc);

   /// scr[s2][] = A(scr[sc1][])
   void multOp(int sc1, int sc2);

   /// scr[sc][] = initial right-hand-side (all subdomains)
   void rhs(int sc);

   void genInv(int e, int type);

   void getCoordNode(int e, int n, ldouble *x);


   void print(const char *s, int sc);

   void print(int sc);

   /// bdValues[][] -= scr[sc][] in all subdomains
   void diffValues(int sc);

   /// bdValues[][] = scr[sc][] from all subdomains
   void fromSubdomains(int sc);

   /// bdValues[][] (primals only) = scr[sc][] (primals)
   void getPrimals(int sc);

   /// scr[sc][] = bdValues  all subdomains
   void setPrimals(int sc);

   /// scr[sc][] = bdValues[][]  all subdomains
   void toSubdomains(int sc);


};

#endif
