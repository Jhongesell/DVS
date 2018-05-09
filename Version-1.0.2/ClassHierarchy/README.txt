Para generar la documentacion del software hay que tener instalado doxygen
en debian para instalar usar:
   # aptitude install doxygen graphviz

Para usarlo hacer:
   $ ./Link
   $ doxygen
   $ ./UnLink

La documentacion HTML quedara en en direcctorio html, para leer la documetacion 
leer el archivo htlm/index.html con cualquier navegador web.



Esta carpeta tiene los siguientes archivos y carpetas

README.txt              Esta documentacion
main.hpp                Archivo descriptivo de la documentacion para doxygen
Doxyfile                Configuracion de la documentacion por doxygen
Link                    Genera las ligas simbolicas a los archivos fuentes
UnLink                  Borra las ligas simbolicas a los archivos fuentes
Examples                Carpeta para poner los ejemplos de las clases que lo ameriten
ActualizaEstiloCodigo   Actualiza el estilo del codigo



Doxygen generara las carpetas
html                    Contiene la salida en formato html
xml                     Contiene la salida en formato xml
latex                   Contiene la salida en formato latex (para compilar y generar PDF usar: make pdf)

