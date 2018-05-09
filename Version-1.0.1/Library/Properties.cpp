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



#include "Properties.hpp"
#include <string.h>
#include <stdio.h>

char *Properties::getProperty(const char *s, const char *val)
{
   ce.nameClassFunct("Properties", "getProperty");

   int i;
   char *r = new (nothrow) char[100];
   if (r == NULL) ce.memoryError("r");
   for (i = 0; i < 100; i++)
   {
      if (keys[i].length() == 0)
      {
         strcpy(r, val);
         return r;
      }
      if (keys[i].compare(s) == 0)
      {
         strcpy(r, vals[i].c_str());
         return r;
      }
   }
   return r;
}


const char *Properties::setProperty(const char *k, const char *v)
{
   int i;
   for (i = 0; i < 100; i++)
   {
      if (keys[i].length() == 0)
      {
         keys[i] = string(k);
         vals[i] = string(v);
         return "";
      }
      if (keys[i].compare(k) == 0)
      {
         string val = vals[i];
         vals[i] = string(v);
         const char *r = val.c_str();
         return r;
      }
   }
   return "";
}


void Properties::list(void)
{
   size_t i;
   const char *r, *s;
   printf("\nList of properties\n");
   for (i = 0; i < keys.size(); i++)
   {
      r = keys[i].c_str();
      s = vals[i].c_str();
      if (*r == '\0') break;
      printf("\n%s %s", r, s);
   }
   printf("\n\n");
}

void Properties::load(istream& stream)
{
   int i, c;
   string k, v;

   try
   {
      // Cycle until the end of the stream
      while ((c = stream.get()) != EOF)
      {
         // ------------- To get a key -------------
         // Cycle to find a character
         while (isspace(c) && (c = stream.get()) != EOF);
         // End of the cycle when the end of the stream is reached
         if (c == EOF)
            break;
         // The string k is cleaned
         k = "";
         // The string k is filled with a word of the stream
         do
         {
            if (c > 32) k += (char)c;
         }
         while (!isspace(c) && (c = stream.get()) != EOF);
         // End of the cycle when the end of the stream is reached
         if (c == EOF)
            break;
         // ------------- To get a value -------------
         while (isspace(c) && (c = stream.get()) != EOF);
         // The string v is cleaned
         v = "";
         // The string v is filled with a word of the stream
         do
         {
            if (c > 32) v += (char)c;
         }
         while (!isspace(c) && (c = stream.get()) != EOF);
         // The values of k and v are asigned to keys and values
         for (i = 0; i < 100; i++)
            // Until find free space
            if (keys[i].length() == 0)
            {
               keys[i] = k;
               vals[i] = v;
               break;
            }
         // End of the cycle when the end of the stream is reached
         if (c == EOF)
            break;
      }
   }
   catch (...)
   {
      cout << "Error in the function load of the class Properties\n";
   }
}
