type src\matrix.cpp^
 src\accel.cpp^
 src\accelHarmonic.cpp^
 src\accelPointMass.cpp^
 src\azelpa.cpp^
 src\cheb3d.cpp^
 src\DEInteg.cpp^
 src\eccAnom.cpp^
 src\eqnEquinox.cpp^
 src\frac.cpp^
 src\g_accelHarmonic.cpp^
 src\gast.cpp^
 src\ghaMatrix.cpp^
 src\global.cpp^
 src\gmst.cpp^
 src\iers.cpp^
 src\position.cpp^
 src\JPL_Eph.cpp^
 src\legendre.cpp^
 src\LTC.cpp^
 src\meanObliquity.cpp^
 src\measUpdate.cpp^
 src\mjday_tdb.cpp^
 src\mjday.cpp^
 src\nutAngles.cpp^
 src\nutMatrix.cpp^
 src\poleMatrix.cpp^
 src\precMatrix.cpp^
 src\r_x.cpp^
 src\r_y.cpp^
 src\r_z.cpp^
 src\sign.cpp^
 src\timediff.cpp^
 src\timeUpdate.cpp^
 src\varEqn.cpp^
 test\tests.cpp^
 > test\result.cpp  2>nul

g++ test/result.cpp -lm -std=c++23 -o bin/tests.exe
del test\result.cpp
cd bin
tests.exe
pause