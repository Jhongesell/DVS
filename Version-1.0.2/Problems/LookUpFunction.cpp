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



#include "LookUpFunction.hpp"
#include "Constant.hpp"
#include "SinPix.hpp"
#include "SinPixSinPiy.hpp"
#include "SinPinxSinPiny.hpp"
#include "SinPinxSinPinySinPinz.hpp"
#include "SinPixCosPiy.hpp"
#include "ExpXY.hpp"
#include "fExpXY.hpp"
#include "NSfExpXY.hpp"
#include "ExpX.hpp"
#include "SinPiXSinPiYSinPiZ.hpp"
#include "ExpVXY.hpp"
#include "ExpVXYZ.hpp"
#include "ExpXYZ.hpp"
#include "NSfExpXYZ.hpp"
#include "SfExpXYZ.hpp"
#include "Disc.hpp"
#include "Disc13.hpp"
#include "Disc14.hpp"
#include "Disc15.hpp"
#include <string.h>
#include <stdio.h>


FunctionV1 *LookUpFunction::getF(char *s)
{
   ce.nameClassFunct("LookUpFunction", "getF");

   printf("\nLookUp %s", s);

   if (strcmp(s, "const") == 0)
   {
      Constant *f = new (nothrow) Constant(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "Disc") == 0)
   {
      Disc *f = new (nothrow) Disc(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "Disc13") == 0)
   {
      Disc13 *f = new (nothrow) Disc13(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "Disc14") == 0)
   {
      Disc14 *f = new (nothrow) Disc14(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "Disc15") == 0)
   {
      Disc15 *f = new (nothrow) Disc15(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPix") == 0)
   {
      SinPix *f = new (nothrow) SinPix(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPixSinPiy") == 0)
   {
      SinPixSinPiy *f = new (nothrow) SinPixSinPiy(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPixCosPiy") == 0)
   {
      SinPixCosPiy *f = new (nothrow) SinPixCosPiy(2.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "ExpXY") == 0)
   {
      ExpXY *f = new (nothrow) ExpXY(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "fExpXY") == 0)
   {
      fExpXY *f = new (nothrow) fExpXY(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "NSfExpXY") == 0)
   {
      NSfExpXY *f = new (nothrow) NSfExpXY(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "ExpX") == 0)
   {
      ExpX *f = new (nothrow) ExpX(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPixSinPiySinPiz") == 0)
   {
      SinPiXSinPiYSinPiZ *f = new (nothrow) SinPiXSinPiYSinPiZ(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "ExpVXY") == 0)
   {
      ExpVXY *f = new (nothrow) ExpVXY(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "ExpVXYZ") == 0)
   {
      ExpVXYZ *f = new (nothrow) ExpVXYZ(1.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "ExpXYZ") == 0)
   {
      ExpXYZ *f = new ExpXYZ(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SfExpXYZ") == 0)
   {
      SfExpXYZ *f = new (nothrow) SfExpXYZ(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "NSfExpXYZ") == 0)
   {
      NSfExpXYZ *f = new (nothrow) NSfExpXYZ(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPinxSinPiny") == 0)
   {
      SinPinxSinPiny *f = new (nothrow) SinPinxSinPiny(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPinxSinPiny") == 0)
   {
      SinPinxSinPiny *f = new (nothrow) SinPinxSinPiny(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }
   else if (strcmp(s, "SinPinxSinPinySinPinz") == 0)
   {
      SinPinxSinPinySinPinz *f = new (nothrow) SinPinxSinPinySinPinz(0.0);
      if (f == NULL) ce.memoryError("f");
      return f;
   }

   FunctionV1 *f = NULL;
   return f;
}

