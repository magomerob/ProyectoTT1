g++ test/tests.cpp src/*.cpp -lm -std=c++23 -o bin/test.exe

cd bin
test.exe
pause