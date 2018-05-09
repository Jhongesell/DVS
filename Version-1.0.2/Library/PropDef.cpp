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


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "Definiciones.hpp"
#include "PropDef.hpp"


int PropDef::parse(string &file)
{
   if (access(file.c_str(), R_OK))
   {
      cout << "The file couldn't be opened\n";
      exit(1);
   }
   try
   {
      fstream datafile; // Stream to manage the file
      // The file is open
      datafile.open(file.c_str(), fstream::in);
      // The stream is sent to load
      load(datafile);
      // The file is closed
      datafile.close();
   }
   catch (...)
   {
      cout << "The file couldn't be opened\n";
      return 1;
   }

   return 0;
}

int PropDef::parse(int nargs, char *args[])
{
   int iarg = 1; // Index for the arguments
   string aux, key, value, name;
   fstream datafile; // Stream to manage the file

   // Cycle through all the arguments
   while (nargs > 1)
   {
      // Case: key is found
      if (args[iarg][0] == '-')
      {
         value = "";
         aux = args[iarg];
         // The key is gotten
         key = aux.substr(1);
         // It's incremented the index and the nargs is reduced
         iarg++;
         nargs--;
         // Case: value is found
         if (nargs > 0 && args[iarg][0] != '-')
         {
            // The value is gotten
            value = args[iarg];
            // It's incremented the index and the nargs is reduced
            iarg++;
            nargs--;
         }
         // The strings found are set
         setProperty(key.c_str(), value.c_str());
      }
      else if (strcmp(args[iarg], "file") == 0)
      {
         // The name of the file is gotten
         name = args[iarg + 1];
         // Print the number of the file
         cout << "file = "  << name << "\n";
         if (access(name.c_str(), R_OK))
         {
            cout << "The file couldn't be opened\n";
            exit(1);
         }
         // It's incremented the index and the nargs is reduced
         iarg += 2;
         nargs -= 2;
         try
         {
            // The file is open
            datafile.open(name.c_str());
            // The stream is sent to load
            load(datafile);
            // The file is closed
            datafile.close();
         }
         catch (...)
         {
            cout << "The file couldn't be opened\n";
            return 1;
         }
      }
      else
         return 2;
   }

   return 0;
}


ldouble PropDef::getDouble(const char *key, ldouble value)
{
   char v[100];
   char tmp[100];
   const char *p = getProperty(key);
   strcpy(v, p);
   delete []p;
   if (*v == '\0')
   {
#ifdef __Double__
      sprintf(tmp, "%f", value);
#else
      sprintf(tmp, "%Lf", value);
#endif
      setProperty(key, tmp);
      return value;
   }
   return (ldouble) atof(v);
}

ldouble PropDef::getDouble(const char *key)
{
   return getDouble(key, 0.0);
}

int PropDef::getInt(const char *key, int value)
{
   char v[100];
   char tmp[100];
   const char *p = getProperty(key);
   strcpy(v, p);
   delete []p;
   if (*v == '\0')
   {
      sprintf(tmp, "%d", value);
      setProperty(key, tmp);
      return value;
   }
   return atoi(v);
}

int PropDef::getInt(const char *key)
{
   return (int) getDouble(key);
}

char *PropDef::getString(const char *key, const char *value)
{
   ce.nameClassFunct("PropDef", "getString");

   char *s = new (nothrow) char[100];
   if (s == NULL) ce.memoryError("s");
   const char *p = getProperty(key);
   strcpy(s, p);
   delete []p;
   if (*s == '\0')
   {
      setProperty(key, value);
      strcpy(s, value);
   }
   return s;
}

const char *PropDef::getString(const char *key)
{
   return getProperty(key);
}

