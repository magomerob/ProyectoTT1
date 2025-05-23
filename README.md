# Proyecto TT1 

Proyecto final de Taller Transversal I.
Este ha consistido en la traducción del codigo en Matlab creado por Meysam Mahooti.

El comando para compilar el archivo principal es:
```
g++ test/EKF_GEOS3.cpp src/*.cpp -lm -std=c++23 -o bin/EKF_GEOS3.exe
```
Este está incluido en el archivo `makefile.bat`, que compila y ejecuta el resultado directamente.

El comando para compilar el fichero de tests es:
```
g++ test/tests.cpp src/*.cpp -lm -std=c++23 -o bin/test.exe
```
Este está incluido en el archivo `makefile_test.bat`, que compila y ejecuta el resultado directamente.