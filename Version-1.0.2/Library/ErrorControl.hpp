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


#ifndef __ErrorControl__
#define __ErrorControl__


#include <new>

using namespace std;


#include <stdlib.h>
#include <stdio.h>


/// Error Control.
/**
 * @author Antonio Carrillo
 * @date Winter 2010
 * @version 0.0.1
 * @bug No errors detected
 * @todo Exception handling
 */
class ErrorControl
{
private:

   /// Name of class
   const char *nmClass;
   /// Name of function
   const char *nmFunction;




public:


   /**
    * Class Constructor
    */
   ErrorControl(void);

   /**
    * Class Constructor
    * @param clas Class name
    */
   ErrorControl(const char *clas);


   /**
    * Class Constructor
    * @param clas Class name
    * @param fun Function name
   */
   ErrorControl(const char *clas, const char *fun);


   /**
    * Name of class and function
    * @param clas Class name
    * @param func Function name
    */
   void nameClassFunct(const char * clas, const char *func);

   /**
    * No memory for this request
    * @param var Var name
    */
   void memoryError(const char * var);

   /**
    * No memory for this request
    * @param var Var name
    * @param i Index number
    */
   void memoryError(const char * var, int i);


   /**
    * No memory for this request
    * @param var  Var name
    * @param func Function name
    */
   void memoryError(const char * var, const char *func);

   /**
    * Fatal error.
    * @param cod Error code
    */
   void fatalError(int cod);

   /**
    * Fatal error.
    * @param cod Error code
    * @param txt Text for user
    */
   void fatalError(int cod, const char *txt);

   /**
    * Set name of class
    * @param clas Class name
    */
   void nameClass(const char *clas);

   /**
    * Set name of function
    * @param func Function name
    */
   void nameFunct(const char *func);

};


/**
 * Error Control.
 *
 * @example ExampleErrorControl.cpp
 */

#endif


