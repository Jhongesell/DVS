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



#ifndef __PropDef__
#define __PropDef__


#include <string>
#include <fstream>
#include <stdarg.h>
#include "Definiciones.hpp"
#include "Properties.hpp"



class PropDef : public Properties
{

public:

   PropDef(void) : Properties()
   {
   }

   PropDef(Properties prop) : Properties(prop)
   {
   }


   PropDef(int nargs, char *args[])  :  Properties()
   {
      // The parse is made
      parse(nargs, args);
   }

   int parse(string &file);

   int parse(int nargs,  char *args[]);

   ldouble getDouble(const char *key, ldouble value);

   ldouble getDouble(const char *key);

   int getInt(const char *key, int value);

   int getInt(const char *key);

   char *getString(const char *key, const char *value);

   const char *getString(const char *key);


};

#endif


