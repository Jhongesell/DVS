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
#include <math.h>
#include "Definiciones.hpp"
#include "MatrizDispersa.hpp"




// Solicita la memoria necesaria para contener los valores de la matriz
// @param ren Numero de renglones de la matriz
// @param col Numero de columnas de la matriz
// @param ban Tamano de la banda
void MatrizDispersa::solicitaMemoria(const int ren, const int col, const int ban)
{
   ce.nameClassFunct("MatrizDispersa", "solicitaMemoria");

   int i;

   // Tamano de la matriz
   Ren = ren;
   Col = col;
   Ban = ban;// > col? col : ban;

   // Solicita memoria para el contenido de la matriz
   M = new (nothrow) ldouble *[Ren];
   if (M == NULL) ce.memoryError("M");
   for (i = 0; i < Ren; i++) M[i] = NULL;
   M[0] = new (nothrow) ldouble[Ban];
   if (M[0] == NULL) ce.memoryError("M[0]", 0);
   for (i = 0; i < Ban; i++) M[0][i] = 0.0;


   // Solicita memoria para almacenar la columna dentro de la banda de la matriz
   J = new (nothrow) int *[Ren];
   if (J == NULL) ce.memoryError("J");
   for (i = 0; i < Ren; i++) J[i] = NULL;
   J[0] = new (nothrow) int[Ban];
   if (J[0] == NULL) ce.memoryError("J[0]", 0);
   for (i = 0; i < Ban; i++) J[0][i] = -1;
}



// Asigna el valor indicado en el renglo y columna solicitado
// @param ren Renglon
// @param col Columna
// @param val Valor
void MatrizDispersa::asigna(const int ren, const int col, const ldouble val)
{
   ce.nameClassFunct("MatrizDispersa", "asigna");

#ifdef DEPURAR
   if (ren < 0 || ren >= Ren || col < 0 || col >= Col)
   {
      printf("\nError al asignar, indices desbordados (%d, %d)\n", ren, col);
      exit(1);
   }
#endif

   int i, k;

   // Solicita la memoria necesaria para almacenar los datos del renglon
   if (M[ren] == NULL)
   {
      M[ren] = new (nothrow) ldouble[Ban];
      if (M[ren] == NULL) ce.memoryError("M[ren]", ren);
      J[ren] = new (nothrow) int[Ban];
      if (J[ren] == NULL) ce.memoryError("J[ren]", ren);
      for (i = 0; i < Ban; i++)
      {
         M[ren][i] = 0.0;
         J[ren][i] = -1;
      }
   }

   // Asignacion de un valor igual a cero dentro de la matriz bandada
   if (val == 0.0)
   {
      // Checar que no existan valores en esa posicion, en caso de existir, eliminar este valor
      k = 0;
      while (k < Ban)
      {
         if (J[ren][k] == -1) break;
         // Columna encontrada
         if (J[ren][k] == col)
         {
            // Se reacomodan los valores de la banda
            for (i = k + 1; i < Ban; i++)
            {
               M[ren][i - 1] = M[ren][i];
               J[ren][i - 1] = J[ren][i];
            }
            M[ren][Ban - 1] = 0.0;
            J[ren][Ban - 1] = -1;
            return;
         }
         // Busca la siguiente columna
         k++;
      }
   }
   else
   {
      // Asignacion de un valor distinto de cero dentro de la matriz bandada
      k = 0;
      while (J[ren][k] != -1 && k < Ban)
      {
         // Columna encontrada
         if (J[ren][k] == col)
         {
            // Se cambia el valor
            M[ren][k] = val;
            return;
         }
         // Busca la siguiente columna
         k++;
#ifdef DEPURAR
         if (k >= Ban)
         {
            printf("\nError al asignar, banda desbordada (%d, %d)\n", ren, col);
            exit(1);
         }
#endif
      }
      // Se introduce el valor en el primer lugar libre encontrado
      M[ren][k] = val;
      J[ren][k] = col;
   }
}



// Retorna el numero de columna cuando se para en el renglon e indice de la banda
// @param ren Numero de renglon
// @param ind Numero de indice
// @return Numero de columna cuando se para en el renglon e indice de la banda
ldouble MatrizDispersa::retorna(const int ren, const int col)
{
#ifdef DEPURAR
   if (ren < 0 || ren >= Ren || col < 0 || col >= Col)
   {
      printf("\nError al recuperar, indices desbordados (%d, %d)\n", ren, col);
      exit(1);
   }
#endif

   if (M[ren] == NULL) return 0.0;

   int k = 0;
   // Busqueda de la columna
   while (J[ren][k] != -1)
   {
      // Columna encontrada
      if (J[ren][k] == col)
      {
         return (M[ren][k]);
      }
      k++;
      if (k >= Ban) break;
   }

   // Columna no encontrada, se asume que el valor guardado es cero
   return 0.0;
}


// Retorna el numero de columnas de la banda para el renglon indicado
// @param ren Numero de renglon
// @return Numero de columnas de la banda para el renglon solicitado
int MatrizDispersa::retornaNumeroColumnasBanda(int ren)
{
#ifdef DEPURAR
   if (ren < 0 || ren >= Ren)
   {
      printf("\nError al buscar el renglon solicitado, indice desbordado %d\n", ren);
      exit(1);
   }
#endif

   if (M[ren] == NULL) return 0;
   int k = 0;
   // Retorna el tamano de la banda
   while (k < Ban)
   {
      if (J[ren][k] == -1) break;
      k++;
   }
   return k;
}





// Multiplica la matriz por el vector B dejando el Resultado en R
// @param b Puntero a un Vector
// @param r Puntero a un Vector
void MatrizDispersa::multiplica(ldouble *b, ldouble *r)
{
   int i, k;
   ldouble v;

   for (i = 0; i < Ren; i++)
   {
      v = 0.0;
      k = 0;
      while (J[i][k] != -1)
      {
         v += M[i][k] *  b[J[i][k]];
         k++;
         if (k >= Ban) break;
      }
      r[i] = v;
   }
}

// Visualiza la matriz
// @param tp (1) Se visualiza el vector de en formato de notacion cientifica, (0) formato notacion de punto flotante
void MatrizDispersa::visualiza(const int tp)
{
   int i, j;
   printf("\n");
   if (tp)
   {
      for (i = 0; i < Ren; i++)
      {
         for (j = 0; j < Col; j++)
         {
#ifdef __Double__
            printf(" %e ", retorna(i, j));
#else
            printf(" %Le ", retorna(i, j));
#endif
         }
         printf("\n");
      }
   }
   else
   {
      for (i = 0; i < Ren; i++)
      {
         for (j = 0; j < Col; j++)
         {
#ifdef __Double__
            printf(" %f ", retorna(i, j));
#else
            printf(" %Lf ", retorna(i, j));
#endif
         }
         printf("\n");
      }
   }
   printf("\n");
}


/// Visualiza las matrices internas usadas para soportar a las matrices bandadas
void MatrizDispersa::visualizaMatricesInternas(void)
{
   int i, j;

   // Visualiza la matriz de datos
   printf("\n\nMatriz de Datos\n");
   for (i = 0; i < renglones(); i++)
   {
      if (!M[i]) continue;
      for (j = 0; j < Ban; j++)
      {
#ifdef __Double__
         printf(" %e ", M[i][j]);
#else
         printf(" %Le ", M[i][j]);
#endif
      }
      printf("\n");
   }

   // Visualiza la matriz de indices
   printf("\n\nMatriz de indices\n");
   for (i = 0; i < renglones(); i++)
   {
      if (!M[i]) continue;
      for (j = 0; j < Ban; j++)
      {
         printf(" %d ", J[i][j]);
      }
      printf("\n");
   }
   printf("\n");
}
