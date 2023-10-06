# Fast Fourier Transform
A C++ Implementation of the Fast Fourier Transform. Written using the radix 2 Cooley and Tukey Aglorithm and no recursive form with zero-padding option.

# Math formulas
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}X_{(N)k}=X_{(N/2)(even)k}&plus;W_{N}^{k}X_{(N/2)(odd)k}" title="{\color{White}X_{(N)k}=X_{(N/2)(even)k}+W_{N}^{k}X_{(N/2)(odd)k}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}X_{(N)N/2&plus;k}=X_{(N/2)(even)k}-W_{N}^{k}X_{(N/2)(odd)k}" title="{\color{White}X_{(N)N/2+k}=X_{(N/2)(even)k}-W_{N}^{k}X_{(N/2)(odd)k}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}k=0,1,...,N/2-1}" title="{\color{White}k=0,1,...,N/2-1}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}W_{N}=e^{\frac{-j2\pi}{N}}}" title="{\color{White}W_{N}=e^{\frac{-j2\pi}{N}}}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}X_{(1)0}=x_{0}}" title="{\color{White}X_{(1)0}=x_{0}}" />

# Compilation 
you need at least c++17

```
g++-13 -std=c++17 fourier.cpp -o fourier -DTIMER -DDIFFERENCE
```
```
./fourier 1024
```

Directives:

-DTIMER (turn on the timer for fft function execution)

-DALLPOINTS (show an array of points as a MatLab form array)

-DDIFFERENCE (the difference of standart deviation (or time if -DTIMER flag on) of recursive and no recursive function)

Example of output:

```
fft (recursive):     2950 microseconds
fft (no recursive):  1639 microseconds
standard deviation (recursive):    [8.4395e-15 + -1.27764e-14i]
standard deviation (no recursive): [1.45703e-16 + 1.48374e-16i]
```
how we can see the non-recursive algorithm works faster around 1,8 times


