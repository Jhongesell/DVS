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


#ifndef __MatrizDispersa__
#define __MatrizDispersa__

#include "Definiciones.hpp"
#include "Matriz_Base.hpp"


/// Clase para el trabajar con matrices dispersas de punto flotante basada en el algoritmo Jagged Diagonal Storage (JDS)
/// El algoritmo esta optimizado para hacer producto matriz vector
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
    @todo Hacer comportamiento para cambiar tamano de banda
    @todo Multiplicacion de matrices
*/
class MatrizDispersa : public Matriz_Base
{
private:

   /// Solicita la memoria necesaria para contener los valores de la matriz
   /** @param ren Numero de renglones de la matriz
       @param col Numero de columnas de la matriz
       @param ban Tamano de la banda  */
   void solicitaMemoria(const int ren, const int col, const int ban);



   /// Puntero a la matriz de datos
   ldouble **M;
   /// Arreglo que contiene los columnas de la matriz
   int **J;


public:

   /// Constructor de la clase
   /** @param ren Numero de renglones de la matriz
       @param col Numero de columnas de la matriz
       @param ban Tamano de la banda   */
   MatrizDispersa(const int ren, const int col, const int ban) : Matriz_Base()
   {
      solicitaMemoria(ren, col, ban);
      asignaNombre("SinNombre");
   }

   /// Constructor de la clase
   /** @param ren Numero de renglones de la matriz
       @param col Numero de columnas de la matriz
       @param ban Tamano de la banda
       @param nmb Nombre de la matriz*/
   MatrizDispersa(const int ren, const int col, const int ban, const char *nmb) : Matriz_Base()
   {
      solicitaMemoria(ren, col, ban);
      asignaNombre(nmb);
   }

   ~MatrizDispersa()
   {
      int i;

      // Libera memoria para el contenido de la matriz
      if (M)
      {
         for (i = 0; i < Ren; i++)
         {
            if (M[i])
            {
               delete []M[i];
               M[i] = NULL;
               delete []J[i];
               J[i] = NULL;
            }
         }
         delete []M;
         M = NULL;
      }

      delete []J;
      J = NULL;
   }


   /// Retorna el tamano de la banda
   /** @return Tamano de la banda */
   inline int tamanoBanda(void)
   {
      return Ban;
   }


   /// Inicializa la matriz al valor indicado
   /** @param val Valor por omision para inicializar la matriz */
   void inicializa(ldouble val)
   {
      printf("\nFuncion no implementada en esta clase\n");
   }


   /// Asigna el valor indicado en el renglo y columna solicitado
   /** @param ren Renglon
       @param col Columna
       @param val Valor */
   void asigna(const int ren, const int col, const ldouble val);

   /// Retorna el numero de columna cuando se para en el renglon e indice de la banda
   /** @param ren Numero de renglon
       @param col Numero de columna
       @return Numero de columna cuando se para en el renglon e indice de la banda  */
   ldouble retorna(const int ren, const int col);

   /// Retorna el numero de columnas de la banda para el renglon indicado
   /** @param ren Numero de renglon
       @return Numero de columnas de la banda para el renglon solicitado  */
   int retornaNumeroColumnasBanda(int ren);

   /// Retorna el numero de columna cuando se para en el renglon e indice de la banda
   /** @param ren Numero de renglon
       @param ind Numero de indice
       @return Numero de columna cuando se para en el renglon e indice de la banda  */
   inline int retornaNumeroColumna(int ren, int ind)
   {
      return (J[ren][ind]);
   }

   /// Retorna el valor de la columna cuando se para en el renglon e indice de la banda
   /** @param ren Numero de renglon
       @param ind Numero de indice
       @return Valor de la columna cuando se para en el renglon e indice de la banda  */
   inline ldouble retornaValorColumna(int ren, int ind)
   {
      return (M[ren][ind]);
   }


   /// Multiplica la matriz por el vector B dejando el Resultado en R
   /** @param b Puntero a un Vector
       @param r Puntero a un Vector */
   void multiplica(ldouble *b, ldouble *r);


   /// Visualiza la matriz
   /** @param tp (1) Se visualiza el vector de en formato de notacion cientifica, (0) formato notacion de punto flotante */
   void visualiza(const int tp);

   /// Visualiza las matrices internas usadas para soportar a las matrices bandadas
   void visualizaMatricesInternas(void);

};


/**
Esta clase implementa los componentes para el trabajar con matrices dispersas de punto flotante
@example EjemploMatrizDispersa.cpp */

#endif
