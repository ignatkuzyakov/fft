# Fast Fourier Transform
A C++ Implementation of the Fast Fourier Transform. Written using the radix 2 Cooley and Tukey Aglorithm with zero-padding option.

# Math formulas
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}X_{(N)k}=X_{(N/2)(even)k}&plus;W_{N}^{k}X_{(N/2)(odd)k}" title="{\color{White}X_{(N)k}=X_{(N/2)(even)k}+W_{N}^{k}X_{(N/2)(odd)k}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}X_{(N)N/2&plus;k}=X_{(N/2)(even)k}-W_{N}^{k}X_{(N/2)(odd)k}" title="{\color{White}X_{(N)N/2+k}=X_{(N/2)(even)k}-W_{N}^{k}X_{(N/2)(odd)k}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}k=0,1,...,N/2-1}" title="{\color{White}k=0,1,...,N/2-1}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}W_{N}=e^{\frac{-j2\pi}{N}}}" title="{\color{White}W_{N}=e^{\frac{-j2\pi}{N}}}" />
<img src="https://latex.codecogs.com/svg.image?\inline&space;\LARGE&space;{\color{White}X_{(1)0}=x_{0}}" title="{\color{White}X_{(1)0}=x_{0}}" />

# Compilation 
you need at least c++14

```
g++-13 -std=c++14 fourier.cpp -o fourier
```
```
./fourier 1024
```

Directives:

-DTIMER (turn on the timer for fft function execution)

-DALLERRS (show an array of errors as a MatLab form array)
