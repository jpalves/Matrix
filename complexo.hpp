#ifndef COMPLEXO_HPP
#define COMPLEXO_HPP

#include <complex>
#include <iomanip>
#include <string>
#include <sstream>
#include <concepts>
#include <cmath>
#include <type_traits>

template <typename T>
class complexo : public std::complex<T> {
public:
    complexo(T r = T{}, T i = T{}) : std::complex<T>(r, i) {}
    complexo(const std::complex<T>& c) : std::complex<T>(c) {}

    // Conversão a partir de complexos com outro tipo escalar U (ex.: math::BigFloat para double)
    template <typename U>
    requires (!std::is_same_v<T, U>) && requires(U u) { static_cast<T>(u); }
    complexo(const std::complex<U>& c)
        : std::complex<T>(static_cast<T>(c.real()), static_cast<T>(c.imag())) {}

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator*(S s, const complexo& z) {
        return complexo(s * z.real(), s * z.imag());
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator*(const complexo& z, S s) {
        return complexo(z.real() * s, z.imag() * s);
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator+(S s, const complexo& z) {
        return complexo(s + z.real(), z.imag());
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator+(const complexo& z, S s) {
        return complexo(z.real() + s, z.imag());
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator-(S s, const complexo& z) {
        return complexo(s - z.real(), -z.imag());
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator-(const complexo& z, S s) {
        return complexo(z.real() - s, z.imag());
    }

    template <typename S>
    requires std::is_arithmetic_v<S>
    friend complexo operator/(const complexo& z, S s) {
        return complexo(z.real() / s, z.imag() / s);
    }
};

template <typename U>
std::string format_scalar(const U& val) {
    if constexpr (requires { val.to_string(); }) {
        if (val == U{0}) return "0";
        return val.to_string();
    } else {
        using std::abs;
        if (abs(val) < 1e-12) return "0";
        if (abs(val - std::round(val)) < 1e-9) {
            return std::to_string(static_cast<long long>(std::round(val)));
        }
        std::ostringstream ss;
        ss << std::defaultfloat << std::setprecision(5) << val;
        return ss.str();
    }
}

template <typename U>
std::ostream &operator<<(std::ostream &stream, const complexo<U>& z) {
    using std::abs;
    U r = z.real();
    U im = z.imag();
    if constexpr (requires { abs(r); }) {
        if (abs(r) < 1e-12) r = U{0};
        if (abs(im) < 1e-12) im = U{0};
    }

    bool r_zero = (r == U{0});
    bool im_zero = (im == U{0});

    std::string s;
    if (r_zero && im_zero) {
        s = "0";
    } else if (im_zero) {
        s = format_scalar(r);
    } else if (r_zero) {
        if (im == U{1}) s = "j";
        else if (im == U{-1}) s = "-j";
        else s = format_scalar(im) + "j";
    } else {
        s = format_scalar(r);
        if (im == U{1}) {
            s += "+j";
        } else if (im == U{-1}) {
            s += "-j";
        } else if (im > U{0}) {
            s += "+" + format_scalar(im) + "j";
        } else {
            s += format_scalar(im) + "j";
        }
    }

    return stream << s;
}

template <typename T>
complexo<T> sqrt(const complexo<T>& z) {
    return std::sqrt(static_cast<const std::complex<T>&>(z));
}

template <typename T>
complexo<double> sqrt(T n) requires std::integral<T> {
    return complexo<double>(std::sqrt(static_cast<complexo<double>>(n)));
}

namespace mat {
template <typename T> struct is_custom_scalar;
template <typename T> struct is_custom_scalar<complexo<T>> : std::true_type {};
}

#endif // COMPLEXO_HPP
