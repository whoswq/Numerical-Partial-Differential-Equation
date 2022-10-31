@echo off
echo compile problem 1
g++ -c My_Matrix.cpp -O3
g++ -c 1_poisson_solver.cpp -O3 
g++ My_Matrix.o 1_poisson_solver.o -o 1_poisson_solver.exe -O3
del *.o
pause