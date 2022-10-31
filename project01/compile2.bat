@echo off
echo compile problem 2
g++ -c My_Matrix.cpp -O3
g++ -c 2_diffusion_CN.cpp -O3 
g++ My_Matrix.o 2_diffusion_CN.o -o 2_diffusion_CN.exe -O3
del *.o
pause