#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>

using complex = std::complex<double>;

class Fourier
{
public:
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
            if(count % 2 == 0)
                even[count / 2] = elem;
            else 
                odd[count / 2] = elem;
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
        N++;
        return N;
    }
};

int main()
{
    size_t points;
    const double step = 0.01;

    std::vector<complex> result;

    for (double i = 0; i <= 2 * M_PI; i += step)
        result.push_back({i, sin(i)});

    points = result.size();
    points = Fourier::next_highest_power_of_2(points);

    size_t n = result.size();
    result.reserve(points);

    for (double i = n; i < points; ++i)
        result.push_back(0);

    std::vector<complex> input = result;

    Fourier::fft(std::begin(result), std::end(result));

    std::cout << '[';

    for (auto &&x : result)
    {
        std::cout << x.real() << " + " << x.imag() << 'i';
        std::cout << "; ";
    }
    std::cout << ']';

    Fourier::ifft(std::begin(result), std::end(result));

    int i = 0;
    for (auto &&x : result)
        x -= input[i++];
}