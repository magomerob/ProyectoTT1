g++ test/EKF_GEOS3.cpp src/*.cpp -lm -o bin/EKF_GEOS3.exe -pg -no-pie

cd bin
EKF_GEOS3.exe
gprof EKF_GEOS3.exe gmon.out > profile.txt
pause