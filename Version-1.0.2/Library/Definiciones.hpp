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




#ifndef __DEFINICIONES_HPP__
#define __DEFINICIONES_HPP__


/// Usa la libreria libgmm++-dev
//#define GMM
#ifdef GMM
   //#define SIN_PRECONDICIONADOR
   #define PRECONDICIONADOR_ILDLTT
#endif



/// Numero maximo de iteraiones en los metodos iterativos
#define NMAXITER 5000                 // Iteraciones en metodos globales
#define NMAXITER_LOCAL 50000          // Iteraciones en metodos locales


/// Tolerancia en los metodos iterativos
#define EPSILON  1e-6                 // Tolerancia en metodos globales
#define EPSILON_LOCAL (EPSILON/1e+4)  // Tolerancia en metodos locales


/// Se toman como iguales dos nodos que difieran en menos que esta EPS_EQUAL
#define EPS_EQUAL 1e-15


/// Con esta opcion visualiza o no el residual de cada iteracion
//#define RESIDUAL


/// Dimension del vector (1) escalar
#define DIM_VECTOR 1


/// Con esta opcion se calcula el numero de condicionamiento en los metodos precondicionados
//#define NUMERO_CONDICIONAMIENTO


/// Activar el modo de depuracion
//#define DEPURAR



/// Definiciones Generales, en caso de no existir definicion generales, solo se consideran coeficientes constantes
//#define COEF_CONSTANTES_ESTABILIZA
//#define COEF_VARAIBLES_ESTABILIZA

/// Definicion de problemas que requieren activar codigo particular para cada problema de ejemplo
//#define NoSimetricoVariable
//#define NoSimetricoRotacional
//#define NoSimetricoNeumann


/// Activacion de las diferentes definiciones para cada problema
#if defined COEF_VARAIBLES_ESTABILIZA
#define ESTABILIZA
#undef  COEFICIENTES_CONSTANTES
#undef  USAR_LOCAL_MET_ITER

#elif defined COEF_CONSTANTES_ESTABILIZA
#define ESTABILIZA
#define COEFICIENTES_CONSTANTES
#undef  USAR_LOCAL_MET_ITER

#else
#define COEFICIENTES_CONSTANTES
#undef  ESTABILIZA
#undef USAR_LOCAL_MET_ITER

#endif




////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                            //
// Si la ecuacion tiene coeficientes constantes entonse se utiliza una discretizacion usando  //
// el metodo de diferencias finitas centradas donde se calculan los coeficientes solo una vez //
//#define COEFICIENTES_CONSTANTES                                                             //
//                                                                                            //
//                                                                                            //
// Con esta opcion se utiliza una discretizacion usando el metodo de diferencias finitas      //
// centradas estabilizadas usando el metodo de la difusion artificial para la ecuacion de     //
// adveccion difusion                                                                         //
//#define ESTABILIZA                                                                          //
//                                                                                            //
//                                                                                            //
// Con esta opcion se utilizan metodos iterativos en lugar de directos para resolver los      //
// sistemas lineales locales, esto genera un ahorro importante de RAM y es ideal para         //
// descomposiciones de dominio finas                                                          //
//#define USAR_LOCAL_MET_ITER                                                                 //
//                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////



/// Activada para trabajar con numeros double en caso contrario trabajar con long double
#define __Double__

#ifdef __Double__
/// Define ldouble como double
typedef double ldouble;
#else
/// Define ldouble como long double
typedef long double ldouble;
#endif



#endif
