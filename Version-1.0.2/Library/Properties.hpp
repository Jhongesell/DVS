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



#ifndef __Properties__
#define __Properties__

#include <string>
#include <vector>
#include <iostream>
#include "ErrorControl.hpp"



using namespace std;

class Properties
{

protected:
   vector<string> keys;
   vector<string> vals;
   /// Control de errores
   ErrorControl ce;


public:

   Properties(void) : keys(100), vals(100)
   { }


   char *getProperty(const char *s, const char *val);

   inline const char *getProperty(const char *s)
   {
      return getProperty(s, "");
   }

   const char *setProperty(const char *k, const char *v);

   void list(void);

   void load(istream& stream);

};

#endif

