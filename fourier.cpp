#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <random>
#include <iterator>
#include <chrono>
#include <vector>
#include <string>

using complex = std::complex<double>;

// timer
#ifdef TIMER
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;
#endif

class Fourier
{
public:
    // recursive fft
    template <typename Iter>
    static void fft(Iter first, Iter last)
    {
        auto N = std::distance(first, last);

        if (N <= 1)
            return;

        std::vector<complex> even(N / 2);
        std::vector<complex> odd(N / 2);    

        int count = 0;

        std::for_each(first, last, [&](auto &&elem)
                      {
            if(count % 2 == 0) even[count / 2] = elem;
            else  odd[count / 2] = elem;
            ++count; });

        fft(std::begin(even), std::end(even));
        fft(std::begin(odd), std::end(odd));

        for (size_t k = 0; k < N / 2; ++k)
        {
            auto t = std::polar(1.0, (-2 * M_PI * k) / N) * odd[k];
            *(first + k) = even[k] + t;
            *(first + k + N / 2) = even[k] - t;
        }
    }

    template <typename Iter>
    static void ifft(Iter first, Iter last)
    {
        auto N = std::distance(first, last);

        std::for_each(first, last, [](auto &&elem)
                      { elem = {elem.real(), elem.imag() * -1}; });

        fft(first, last);

        std::for_each(first, last, [&](auto &&elem)
                      { elem = {elem.real(), elem.imag() * -1};
                        elem /= N; });
    }

    static size_t next_highest_power_of_2(size_t N)
    {
        N--;
        N |= N >> 1;
        N |= N >> 2;
        N |= N >> 4;
        N |= N >> 8;
        N |= N >> 16;
        N |= N >> 32;
        N++;
        return N;
    }
};

int main(int argc, char const *argv[])
{
    size_t points;
    if (argc == 1)
        points = 256;
    else
        points = std::stoi(argv[1]);

    std::vector<complex> result(points);

    std::random_device rd;
    std::default_random_engine reng(rd());
    std::uniform_int_distribution<int> dist(0, points);

    std::generate(result.begin(), result.end(), [&dist, &reng]
                  { return complex{(double)dist(reng), sin((double)dist(reng))}; });

#ifdef TIMER
    auto tstart = high_resolution_clock::now();
#endif
    points = Fourier::next_highest_power_of_2(points);

    size_t n = result.size();
    result.reserve(points);

    for (double i = n; i < points; ++i)
        result.emplace_back(0);

    std::vector<complex> input = result;

    Fourier::fft(std::begin(result), std::end(result));
#ifdef TIMER
    auto tfin = high_resolution_clock::now();

    std::cout << "recursive fft:  " << duration_cast<microseconds>(tfin - tstart).count()
              << " microseconds" << std::endl;
#endif

    Fourier::ifft(std::begin(result), std::end(result));

#ifndef ALLERRS
    int i = 0;
    complex standard_deviation{};
    for (auto &&x : result)
    {
        x -= input[i++];
        x = std::pow(x, 2);
        standard_deviation += x;
    }
    standard_deviation /= result.size();
    standard_deviation = std::sqrt(standard_deviation);

    std::cout << "standard deviation: " << '[' << standard_deviation.real() << " + "
              << standard_deviation.imag() << 'i' << ']' << std::endl;

#else
    int i = 0;
    for (auto &&x : result)
        x -= input[i++];

    std::cout << '[';
    for (auto &&x : result)
    {
        std::cout << x.real() << " + " << x.imag() << 'i';
        std::cout << "; ";
    }
    std::cout << ']';
#endif
}