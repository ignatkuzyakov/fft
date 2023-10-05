#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <random>
#include <iterator>
#include <chrono>
#include <vector>
#include <string>
#include <functional>

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
    static void recFFT(Iter first, Iter last)
    {
        auto N = std::distance(first, last);

        if (N <= 1)
            return;

        std::vector<complex> even(N / 2);
        std::vector<complex> odd(N / 2);

        size_t count = 0;

        std::for_each(first, last, [&](auto &&elem)
                      {
            if(count % 2 == 0) even[count / 2] = elem;
            else               odd[count / 2] = elem;
            ++count; });

        recFFT(std::begin(even), std::end(even));
        recFFT(std::begin(odd), std::end(odd));

        for (size_t k = 0; k < N / 2; ++k)
        {
            auto t = std::polar(1.0, (-2 * M_PI * k) / N) * odd[k];
            *(first + k) = even[k] + t;
            *(first + k + N / 2) = even[k] - t;
        }
    }

    // no recursive fft
    template <typename Iter>
    static void norecFFT(Iter first, Iter last)
    {
        auto N = std::distance(first, last);
        size_t lg = std::log2(N);

        for (size_t i = 0, rev = 0; i < N; ++i, rev = 0)
        {
            for (size_t j = 0; j < lg; ++j)
                if (i & (1 << j))
                    rev |= 1 << (lg - 1 - j);

            if (i < rev)
                std::swap(*(first + i), *(first + rev));
        }

        for (size_t sz = 2; sz <= N; sz <<= 1)
        {
            for (size_t i = 0; i < N; i += sz)
            {
                complex w(1);
                auto e = std::polar(1.0, (-2 * M_PI) / sz);
                for (size_t j = 0; j < sz / 2; ++j)
                {
                    auto t = w * *(first + i + j + sz / 2);
                    auto even = *(first + i + j);
                    *(first + i + j) = even + t;
                    *(first + i + j + sz / 2) = even - t;
                    w *= e;
                }
            }
        }
    }

    template <typename Iter, typename FFT>
    static void ifft(Iter first, Iter last, FFT&& fft)
    {
        auto N = std::distance(first, last);

        std::for_each(first, last, [](auto &&elem) 
                    { elem = {elem.real(), elem.imag() * -1}; });

        ::std::invoke(std::forward<FFT>(fft), first, last);

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

namespace utils
{
    template<typename T>
    complex standard_deviation(const T& input, const T& output) 
    {
        size_t i = 0;
        T result = output;
        complex standard_deviation{};
        for (auto &&x : result)
        {
            x -= input[i++];
            x = std::pow(x, 2);
            standard_deviation += x;
        }

        standard_deviation /= result.size();
        standard_deviation = std::sqrt(standard_deviation);

        return standard_deviation;
    }
    
    template<typename T>
    void print_difference(const T& input, const T& output)
    {
        size_t i = 0;
        T result = output;
        for (auto &&x : result)
            x -= input[i++];

        std::cout << '[';
        for (auto &&x : result)
        {
            std::cout << x.real() << " + " << x.imag() << 'i';
            std::cout << "; ";
        }
        std::cout << ']';
    }
}

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

    points = Fourier::next_highest_power_of_2(points);

    size_t n = result.size();
    result.reserve(points);

    for (double i = n; i < points; ++i)
        result.emplace_back(0);

    std::vector<complex> input = result;

#ifdef DIFFERENCE

    std::vector<complex> resultNoRec = result;

    #ifdef TIMER

        auto tstartNoRec = high_resolution_clock::now();

    #endif

    Fourier::norecFFT(std::begin(resultNoRec), std::end(resultNoRec));

    #ifdef TIMER

        auto tfinNoRec = high_resolution_clock::now();

    #endif

#endif

#ifdef TIMER

    auto tstart = high_resolution_clock::now();

#endif

    Fourier::recFFT(std::begin(result), std::end(result));

#ifdef TIMER

    auto tfin = high_resolution_clock::now();

    std::cout << "fft (recursive):     " << duration_cast<microseconds>(tfin - tstart).count()
              << " microseconds" << std::endl;

    #ifdef DIFFERENCE

        std::cout << "fft (no recursive):  " << duration_cast<microseconds>(tfinNoRec - tstartNoRec).count()
                << " microseconds" << std::endl;

    #endif

#endif

    Fourier::ifft(std::begin(result), std::end(result), &Fourier::recFFT<std::vector<complex>::iterator>);

#ifndef ALLERRS

    auto standard_deviation = utils::standard_deviation(input, result);

    std::cout << "standard deviation (recursive):    " << '[' << standard_deviation.real() << " + "
              << standard_deviation.imag() << 'i' << ']' << std::endl;

    #ifdef DIFFERENCE

        Fourier::ifft(std::begin(resultNoRec), std::end(resultNoRec),  &Fourier::norecFFT<std::vector<complex>::iterator>);

        auto standard_deviation_NoRec = utils::standard_deviation(input, resultNoRec);

        std::cout << "standard deviation (no recursive): " << '[' << standard_deviation_NoRec.real() << " + "
              << standard_deviation_NoRec.imag() << 'i' << ']' << std::endl;

#endif
        
#else

    utils::print_difference(input, result);
    
    #ifdef DIFFERENCE

        utils::print_difference(input, resultNoRec);
       
    #endif

#endif

}