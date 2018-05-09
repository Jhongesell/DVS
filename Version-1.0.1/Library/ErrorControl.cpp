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



#include "ErrorControl.hpp"


//
// Class Constructor
//
ErrorControl::ErrorControl(void)
{
   nameClassFunct(" ", " ");
}

//
// Class Constructor
// @param clas Class name
//
ErrorControl::ErrorControl(const char *clas)
{
   nameClassFunct(clas, " ");
}


//
// Class Constructor
// @param clas Class name
// @param fun Function name
//
ErrorControl::ErrorControl(const char *clas, const char *fun)
{
   nameClassFunct(clas, fun);
}


//
// Name of class and function
// @param clas Class name
// @param func Function name
//
void ErrorControl::nameClassFunct(const char * clas, const char *func)
{
   nameClass(clas);
   nameFunct(func);
}

//
// No memory for this request
// @param var Var name
//
void ErrorControl::memoryError(const char * var)
{
   printf("\n\nNo memory for %s request in %s of class %s\n\n", var, nmFunction, nmClass);
   fatalError(1);
}

//
// No memory for this request
// @param var Var name
// @param i Index number
//
void ErrorControl::memoryError(const char * var, int i)
{
   printf("\n\nNo memory for %s request %d in %s of class %s\n\n", var, i, nmFunction, nmClass);
   fatalError(1);
}


//
// No memory for this request
// @param var  Var name
// @param func Function name
//
void ErrorControl::memoryError(const char * var, const char *func)
{
   printf("\n\nNo memory for %s request in %s of class %s\n\n", var, func, nmClass);
   fatalError(1);
}

//
// Fatal error.
// @param cod Error code
//
void ErrorControl::fatalError(int cod)
{
   printf("\nEnd program\n");
   exit(cod);
}

//
// Fatal error.
// @param cod Error code
// @param txt Text for user
//
void ErrorControl::fatalError(int cod, const char *txt)
{
   printf("\n %s\nEnd program\n", txt);
   exit(cod);
}


//
// Set name of class
// @param clas Class name
//
void ErrorControl::nameClass(const char *clas)
{
   nmClass = clas;
}

//
// Set name of function
// @param func Function name
//
void ErrorControl::nameFunct(const char *func)
{
   nmFunction = func;
}
