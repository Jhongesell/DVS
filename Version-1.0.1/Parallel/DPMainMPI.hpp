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



#ifndef __DPMainMPI__
#define __DPMainMPI__


#include <time.h>
#include "Definiciones.hpp"
#include "EsquemaMEMPI.hpp"
#include "PropDef.hpp"
#include "EllipOp.hpp"
#include "InternalBd.hpp"
#include "FunctionV1.hpp"
#include "Primal.hpp"
#include "EllipOp.hpp"
#include "Constant.hpp"
#include "LookUpFunction.hpp"
#include "VertPrimal.hpp"
#include "VertEdgePrimal.hpp"
#include "AllPrimal.hpp"
#include "NoPrimal.hpp"
#include "RectSub.hpp"



/// Clase base para definir a los metodos DVS-DDM
/** Clase base para definir a los metodos DVS-DDM en paralelo
en donde se definen las operaciones que realizaran los nodos
esclavos del esquema Mestro-Esclavo y la inicializacion de la
parte paralela de la ejecucion
*/
/** @author Antonio Carrillo Ledesma
    @date primavera 2010
    @version 1.0.0
    @bug No hay errores conocidos
*/
class DPMainMPI:  public EsquemaMEMPI
{
protected:

   /// Tiempo inicial
   time_t t1;
   /// Tiempo final
   time_t t2;
   /// Nmero de tareas por nodo esclavo
   int nta;
   /// Nmero de esclavo en el que estara la tarea
   int xnp;
   /// Nmero de tarea dentro del esclavo
   int indl;
   /// Arreglo para recibir mensajes
   int msa[10];
   /// Arreglo para enviar mensajes
   int mss[10];



   vector<InternalBd*>  hbd;
   vector<RectSub*>   omegas;



   FunctionV1 *zero, *one;
   FunctionV1 *f, *g;
   char *sf, *sg;
   ldouble fc, gc;
   int *mesh;
   char *prim;
   char *method;
   int swprint;

   ldouble **domain;
   ldouble c;
   Primal *primal;
   EllipOp *op;

   int nDim;
   int nOmega;



   void deleteInternalBd(void);



public:

   /// Constructor de la clase
   DPMainMPI(int id, int np, PropDef &props, EllipOp &op);

   /// Destructor de la clase
   ~DPMainMPI();

   /// Esclavo
   void Esclavo(void);

};


#endif


