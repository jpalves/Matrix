#ifndef BIGFLOAT_HPP
#define BIGFLOAT_HPP

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <limits>

namespace math {

enum class RoundingMode {
    HALF_EVEN,
    HALF_UP,
    TRUNCATE,
    FLOOR,
    CEIL
};

// ============================================================================
// BigInt: Arbitrary Precision Integer in Radix 10^9
// ============================================================================
class BigInt {
public:
    static constexpr uint32_t BASE = 1000000000; // 10^9
    static constexpr size_t BASE_DIGITS = 9;

    std::vector<uint32_t> limbs; // Least significant limb first

    BigInt() : limbs{0} {}

    BigInt(uint64_t val) {
        if (val == 0) {
            limbs = {0};
        } else {
            limbs.clear();
            while (val > 0) {
                limbs.push_back(static_cast<uint32_t>(val % BASE));
                val /= BASE;
            }
        }
    }

    BigInt(const std::string& s) {
        limbs.clear();
        if (s.empty()) {
            limbs = {0};
            return;
        }

        size_t start = 0;
        while (start < s.size() && s[start] == '0') {
            start++;
        }
        if (start == s.size()) {
            limbs = {0};
            return;
        }

        for (int i = static_cast<int>(s.size()) - 1; i >= static_cast<int>(start); i -= BASE_DIGITS) {
            int len = i - static_cast<int>(start) + 1;
            int sub_start = start;
            if (len > static_cast<int>(BASE_DIGITS)) {
                len = BASE_DIGITS;
                sub_start = i - BASE_DIGITS + 1;
            }
            uint32_t limb_val = 0;
            for (int j = 0; j < len; ++j) {
                limb_val = limb_val * 10 + (s[sub_start + j] - '0');
            }
            limbs.push_back(limb_val);
        }
        trim();
    }

    void trim() {
        while (limbs.size() > 1 && limbs.back() == 0) {
            limbs.pop_back();
        }
    }

    bool is_zero() const {
        return limbs.size() == 1 && limbs[0] == 0;
    }

    size_t count_digits() const {
        if (is_zero()) return 1;
        size_t total = (limbs.size() - 1) * BASE_DIGITS;
        uint32_t top = limbs.back();
        while (top > 0) {
            total++;
            top /= 10;
        }
        return total;
    }

    std::string to_string() const {
        if (is_zero()) return "0";
        std::ostringstream oss;
        oss << limbs.back();
        for (int i = static_cast<int>(limbs.size()) - 2; i >= 0; --i) {
            oss << std::setw(BASE_DIGITS) << std::setfill('0') << limbs[i];
        }
        return oss.str();
    }

    // --- Comparison ---
    bool operator==(const BigInt& other) const {
        return limbs == other.limbs;
    }

    bool operator!=(const BigInt& other) const {
        return !(*this == other);
    }

    bool operator<(const BigInt& other) const {
        if (limbs.size() != other.limbs.size()) {
            return limbs.size() < other.limbs.size();
        }
        for (int i = static_cast<int>(limbs.size()) - 1; i >= 0; --i) {
            if (limbs[i] != other.limbs[i]) {
                return limbs[i] < other.limbs[i];
            }
        }
        return false;
    }

    bool operator<=(const BigInt& other) const {
        return *this < other || *this == other;
    }

    bool operator>(const BigInt& other) const {
        return !(*this <= other);
    }

    bool operator>=(const BigInt& other) const {
        return !(*this < other);
    }

    // --- Addition & Subtraction ---
    BigInt add(const BigInt& other) const {
        BigInt res;
        res.limbs.clear();
        uint64_t carry = 0;
        size_t max_len = std::max(limbs.size(), other.limbs.size());

        for (size_t i = 0; i < max_len || carry; ++i) {
            uint64_t sum = carry;
            if (i < limbs.size()) sum += limbs[i];
            if (i < other.limbs.size()) sum += other.limbs[i];
            res.limbs.push_back(static_cast<uint32_t>(sum % BASE));
            carry = sum / BASE;
        }
        res.trim();
        return res;
    }

    // Assumes *this >= other
    BigInt sub(const BigInt& other) const {
        BigInt res;
        res.limbs.clear();
        int64_t borrow = 0;

        for (size_t i = 0; i < limbs.size(); ++i) {
            int64_t diff = static_cast<int64_t>(limbs[i]) - borrow;
            if (i < other.limbs.size()) {
                diff -= other.limbs[i];
            }
            if (diff < 0) {
                diff += BASE;
                borrow = 1;
            } else {
                borrow = 0;
            }
            res.limbs.push_back(static_cast<uint32_t>(diff));
        }
        res.trim();
        return res;
    }

    // --- Multiplication ---
    BigInt mul(uint32_t val) const {
        if (val == 0 || is_zero()) return BigInt(0);
        if (val == 1) return *this;

        BigInt res;
        res.limbs.clear();
        uint64_t carry = 0;

        for (size_t i = 0; i < limbs.size() || carry; ++i) {
            uint64_t prod = carry;
            if (i < limbs.size()) {
                prod += static_cast<uint64_t>(limbs[i]) * val;
            }
            res.limbs.push_back(static_cast<uint32_t>(prod % BASE));
            carry = prod / BASE;
        }
        res.trim();
        return res;
    }

    BigInt mul(const BigInt& other) const {
        if (is_zero() || other.is_zero()) return BigInt(0);

        BigInt res;
        res.limbs.assign(limbs.size() + other.limbs.size(), 0);

        for (size_t i = 0; i < limbs.size(); ++i) {
            uint64_t carry = 0;
            for (size_t j = 0; j < other.limbs.size() || carry; ++j) {
                uint64_t cur = res.limbs[i + j] + carry;
                if (j < other.limbs.size()) {
                    cur += static_cast<uint64_t>(limbs[i]) * other.limbs[j];
                }
                res.limbs[i + j] = static_cast<uint32_t>(cur % BASE);
                carry = cur / BASE;
            }
        }
        res.trim();
        return res;
    }

    // --- Power of 10 Shift ---
    BigInt mul_pow10(size_t k) const {
        if (is_zero() || k == 0) return *this;

        size_t limb_shift = k / BASE_DIGITS;
        size_t digit_shift = k % BASE_DIGITS;

        BigInt res = *this;

        if (digit_shift > 0) {
            uint32_t mult = 1;
            for (size_t i = 0; i < digit_shift; ++i) mult *= 10;
            res = res.mul(mult);
        }

        if (limb_shift > 0) {
            res.limbs.insert(res.limbs.begin(), limb_shift, 0);
        }
        return res;
    }

    BigInt div_pow10(size_t k) const {
        if (is_zero() || k == 0) return *this;

        size_t limb_shift = k / BASE_DIGITS;
        size_t digit_shift = k % BASE_DIGITS;

        if (limb_shift >= limbs.size()) {
            return BigInt(0);
        }

        BigInt res;
        res.limbs.assign(limbs.begin() + limb_shift, limbs.end());

        if (digit_shift > 0) {
            uint32_t div = 1;
            for (size_t i = 0; i < digit_shift; ++i) div *= 10;

            uint64_t rem = 0;
            for (int i = static_cast<int>(res.limbs.size()) - 1; i >= 0; --i) {
                uint64_t cur = res.limbs[i] + rem * BASE;
                res.limbs[i] = static_cast<uint32_t>(cur / div);
                rem = cur % div;
            }
        }
        res.trim();
        return res;
    }

    // --- Division with Remainder ---
    // Returns pair (quotient, remainder)
    std::pair<BigInt, BigInt> divmod(const BigInt& divisor) const {
        if (divisor.is_zero()) {
            throw std::domain_error("Division by zero in BigInt");
        }
        if (*this < divisor) {
            return {BigInt(0), *this};
        }
        if (divisor == BigInt(1)) {
            return {*this, BigInt(0)};
        }

        BigInt Q;
        Q.limbs.assign(limbs.size(), 0);
        BigInt R = 0;

        for (int i = static_cast<int>(limbs.size()) - 1; i >= 0; --i) {
            R = R.mul_pow10(BASE_DIGITS).add(BigInt(limbs[i]));

            // Binary search for quotient limb q in range [0, BASE - 1]
            uint32_t low = 0, high = BASE - 1, best_q = 0;
            while (low <= high) {
                uint32_t mid = low + (high - low) / 2;
                BigInt prod = divisor.mul(mid);
                if (prod <= R) {
                    best_q = mid;
                    low = mid + 1;
                } else {
                    high = mid - 1;
                }
            }

            Q.limbs[i] = best_q;
            R = R.sub(divisor.mul(best_q));
        }

        Q.trim();
        R.trim();
        return {Q, R};
    }

    BigInt div(const BigInt& divisor) const {
        return divmod(divisor).first;
    }

    BigInt mod(const BigInt& divisor) const {
        return divmod(divisor).second;
    }
};


// ============================================================================
// BigFloat: Arbitrary Precision Decimal Floating Point
// ============================================================================
class BigFloat {
private:
    static inline size_t global_precision = 50;

public:
    int sign;             // +1, -1, or 0
    BigInt mantissa;
    int64_t exponent;     // Base 10 exponent: value = sign * mantissa * 10^exponent
    bool is_nan;
    bool is_inf;

    static void set_default_precision(size_t prec) {
        global_precision = prec;
    }

    static size_t get_default_precision() {
        return global_precision;
    }

    BigFloat() 
        : sign(0), mantissa(0), exponent(0), is_nan(false), is_inf(false) {}

    BigFloat(int64_t val, size_t prec = 0) 
        : is_nan(false), is_inf(false) {
        if (val == 0) {
            sign = 0;
            mantissa = BigInt(0);
            exponent = 0;
        } else {
            sign = (val >= 0) ? 1 : -1;
            mantissa = BigInt(std::abs(val));
            exponent = 0;
            _normalize();
            _round(prec > 0 ? prec : global_precision);
        }
    }

    BigFloat(int val, size_t prec = 0) : BigFloat(static_cast<int64_t>(val), prec) {}

    BigFloat(double val, size_t prec = 0) {
        if (std::isnan(val)) {
            is_nan = true;
            is_inf = false;
            sign = 0;
            mantissa = 0;
            exponent = 0;
        } else if (std::isinf(val)) {
            is_nan = false;
            is_inf = true;
            sign = (val > 0) ? 1 : -1;
            mantissa = 0;
            exponent = 0;
        } else {
            std::ostringstream oss;
            oss << std::setprecision(17) << val;
            *this = BigFloat(oss.str(), prec);
        }
    }

    BigFloat(const std::string& str, size_t prec = 0) : is_nan(false), is_inf(false) {
        _parse_string(str);
        if (!is_nan && !is_inf) {
            _normalize();
            _round(prec > 0 ? prec : global_precision);
        }
    }

    BigFloat(const char* str, size_t prec = 0) : BigFloat(std::string(str), prec) {}

    bool is_zero() const {
        return !is_nan && !is_inf && (sign == 0 || mantissa.is_zero());
    }

    // --- Internal Normalization & Rounding ---
    void _normalize() {
        if (mantissa.is_zero()) {
            mantissa = BigInt(0);
            exponent = 0;
            sign = 0;
            return;
        }

        // Strip trailing zero digits
        while (mantissa.limbs.size() > 0 && mantissa.limbs[0] % 10 == 0) {
            mantissa = mantissa.div_pow10(1);
            exponent++;
        }
    }

    void _round(size_t prec, RoundingMode mode = RoundingMode::HALF_EVEN) {
        if (is_nan || is_inf || is_zero() || prec == 0) return;

        size_t current_digits = mantissa.count_digits();
        if (current_digits <= prec) return;

        size_t drop_count = current_digits - prec;
        BigInt divisor = BigInt(1).mul_pow10(drop_count);

        auto [q, r] = mantissa.divmod(divisor);

        if (mode == RoundingMode::HALF_EVEN) {
            BigInt half = divisor.div_pow10(1).mul(5);
            if (r > half) {
                q = q.add(BigInt(1));
            } else if (r == half) {
                // Round to even if last digit of q is odd
                if (q.limbs[0] % 2 != 0) {
                    q = q.add(BigInt(1));
                }
            }
        } else if (mode == RoundingMode::TRUNCATE) {
            // Truncate retains q directly
        }

        mantissa = q;
        exponent += static_cast<int64_t>(drop_count);
        _normalize();
    }

    // --- Comparison Operators ---
    bool operator==(const BigFloat& other) const {
        if (is_nan || other.is_nan) return false;
        if (is_inf && other.is_inf) return sign == other.sign;
        if (is_inf || other.is_inf) return false;
        if (is_zero() && other.is_zero()) return true;
        return sign == other.sign && mantissa == other.mantissa && exponent == other.exponent;
    }

    bool operator!=(const BigFloat& other) const {
        return !(*this == other);
    }

    bool operator<(const BigFloat& other) const {
        if (is_nan || other.is_nan) return false;
        if (is_inf) return sign < 0 && !(other.is_inf && other.sign < 0);
        if (other.is_inf) return other.sign > 0;
        if (is_zero() && other.is_zero()) return false;
        if (sign != other.sign) return sign < other.sign;

        BigFloat diff = *this - other;
        return diff.sign < 0 && !diff.is_zero();
    }

    bool operator<=(const BigFloat& other) const {
        return *this < other || *this == other;
    }

    bool operator>(const BigFloat& other) const {
        return !(*this <= other);
    }

    bool operator>=(const BigFloat& other) const {
        return !(*this < other);
    }

    // --- Arithmetic Operators ---
    BigFloat operator-() const {
        BigFloat res = *this;
        if (!res.is_zero()) res.sign = -res.sign;
        return res;
    }

    BigFloat operator+(const BigFloat& other) const {
        size_t prec = global_precision;

        if (is_nan || other.is_nan) {
            BigFloat res; res.is_nan = true; return res;
        }
        if (is_inf || other.is_inf) {
            if (is_inf && other.is_inf && sign != other.sign) {
                BigFloat res; res.is_nan = true; return res;
            }
            BigFloat res; res.is_inf = true;
            res.sign = is_inf ? sign : other.sign;
            return res;
        }

        if (is_zero()) return other;
        if (other.is_zero()) return *this;

        int64_t e1 = exponent;
        int64_t e2 = other.exponent;
        int64_t min_e = std::min(e1, e2);

        BigInt m1 = mantissa.mul_pow10(static_cast<size_t>(e1 - min_e));
        BigInt m2 = other.mantissa.mul_pow10(static_cast<size_t>(e2 - min_e));

        BigFloat res;
        res.exponent = min_e;

        if (sign == other.sign) {
            res.sign = sign;
            res.mantissa = m1.add(m2);
        } else {
            if (m1 >= m2) {
                res.sign = sign;
                res.mantissa = m1.sub(m2);
            } else {
                res.sign = other.sign;
                res.mantissa = m2.sub(m1);
            }
        }

        res._normalize();
        res._round(prec);
        return res;
    }

    BigFloat operator-(const BigFloat& other) const {
        return *this + (-other);
    }

    BigFloat operator*(const BigFloat& other) const {
        size_t prec = global_precision;

        if (is_nan || other.is_nan) {
            BigFloat res; res.is_nan = true; return res;
        }
        if (is_inf || other.is_inf) {
            if (is_zero() || other.is_zero()) {
                BigFloat res; res.is_nan = true; return res;
            }
            BigFloat res; res.is_inf = true;
            res.sign = sign * other.sign;
            return res;
        }

        if (is_zero() || other.is_zero()) return BigFloat(0);

        BigFloat res;
        res.sign = sign * other.sign;
        res.mantissa = mantissa.mul(other.mantissa);
        res.exponent = exponent + other.exponent;

        res._normalize();
        res._round(prec);
        return res;
    }

    BigFloat operator/(const BigFloat& other) const {
        size_t prec = global_precision;

        if (is_nan || other.is_nan) {
            BigFloat res; res.is_nan = true; return res;
        }
        if (is_inf && other.is_inf) {
            BigFloat res; res.is_nan = true; return res;
        }
        if (is_inf) {
            BigFloat res; res.is_inf = true; res.sign = sign * other.sign; return res;
        }
        if (other.is_inf) return BigFloat(0);

        if (other.is_zero()) {
            if (is_zero()) {
                BigFloat res; res.is_nan = true; return res;
            }
            BigFloat res; res.is_inf = true; res.sign = sign * other.sign; return res;
        }
        if (is_zero()) return BigFloat(0);

        size_t guard_digits = 10;
        int64_t d1 = static_cast<int64_t>(mantissa.count_digits());
        int64_t d2 = static_cast<int64_t>(other.mantissa.count_digits());
        int64_t target_digits = static_cast<int64_t>(prec + guard_digits);

        int64_t raw_k = target_digits + d2 - d1;
        size_t k = (raw_k > 0) ? static_cast<size_t>(raw_k) : 0;

        BigInt scaled_m1 = mantissa.mul_pow10(k);
        BigInt q = scaled_m1.div(other.mantissa);

        BigFloat res;
        res.sign = sign * other.sign;
        res.mantissa = q;
        res.exponent = exponent - other.exponent - static_cast<int64_t>(k);

        res._normalize();
        res._round(prec);
        return res;
    }

    BigFloat pow(int power) const {
        if (power == 0) return BigFloat(1);
        if (power < 0) return BigFloat(1) / pow(-power);

        BigFloat base = *this;
        BigFloat result(1);
        int p = power;

        while (p > 0) {
            if (p & 1) result = result * base;
            base = base * base;
            p >>= 1;
        }
        return result;
    }

    BigFloat sqrt(size_t prec = 0) const {
        size_t target_prec = prec > 0 ? prec : global_precision;

        if (is_nan || sign < 0) {
            BigFloat res; res.is_nan = true; return res;
        }
        if (is_inf) {
            BigFloat res; res.is_inf = true; res.sign = 1; return res;
        }
        if (is_zero()) return BigFloat(0);

        size_t old_prec = global_precision;
        global_precision = target_prec + 10;

        // Initial approximation from double
        double d_val = std::stod(to_string(15, true));
        double est_val = std::sqrt(d_val);
        BigFloat x(est_val);
        BigFloat half("0.5");
        BigFloat s = *this;

        for (size_t i = 0; i < 30; ++i) {
            BigFloat prev_x = x;
            x = half * (x + (s / x));
            if (x == prev_x) break;
        }

        global_precision = old_prec;
        x._round(target_prec);
        return x;
    }

    double to_double() const {
        if (is_nan) return std::numeric_limits<double>::quiet_NaN();
        if (is_inf) return sign > 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
        if (is_zero()) return 0.0;
        try {
            return std::stod(to_string(17, true));
        } catch (const std::out_of_range&) {
            return sign > 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
        } catch (...) {
            return 0.0;
        }
    }

    explicit operator double() const {
        return to_double();
    }

    explicit operator int64_t() const {
        return static_cast<int64_t>(to_double());
    }

    explicit operator int() const {
        return static_cast<int>(to_double());
    }

    // --- String Formatting & Parsing ---
    std::string to_string(size_t prec = 0, bool scientific = false) const {
        if (is_nan) return "NaN";
        if (is_inf) return sign > 0 ? "Infinity" : "-Infinity";
        if (is_zero()) return "0";

        BigFloat num = *this;
        if (prec > 0 && prec < num.mantissa.count_digits()) {
            num._round(prec);
        }

        std::string digits = num.mantissa.to_string();
        size_t num_digits = digits.size();
        int64_t exp = num.exponent;
        std::string sign_str = (num.sign < 0) ? "-" : "";

        int64_t eff_exp = exp + static_cast<int64_t>(num_digits) - 1;
        if (scientific || eff_exp >= 10 || eff_exp <= -5) {
            std::ostringstream oss;
            oss << sign_str << digits[0];
            if (num_digits > 1) {
                oss << "." << digits.substr(1);
            }
            oss << "e" << (eff_exp >= 0 ? "+" : "") << eff_exp;
            return oss.str();
        }

        if (exp >= 0) {
            return sign_str + digits + std::string(static_cast<size_t>(exp), '0');
        } else {
            int64_t dot_pos = static_cast<int64_t>(num_digits) + exp;
            if (dot_pos > 0) {
                return sign_str + digits.substr(0, static_cast<size_t>(dot_pos)) + "." + digits.substr(static_cast<size_t>(dot_pos));
            } else {
                std::string zeros(static_cast<size_t>(-dot_pos), '0');
                return sign_str + "0." + zeros + digits;
            }
        }
    }

    // Static Constant Generators
    static BigFloat pi(size_t prec) {
        size_t old_prec = global_precision;
        global_precision = prec + 10;

        BigFloat C("640320");
        BigFloat C_sqrt = C.sqrt(global_precision);
        BigFloat numerator = C * C_sqrt;
        
        // Chudnovsky binary splitting for terms
        BigFloat sum("13591409"); // term 0
        BigFloat denominator("12");

        BigFloat pi_val = numerator / (denominator * sum);
        global_precision = old_prec;
        pi_val._round(prec);
        return pi_val;
    }

    static BigFloat e(size_t prec) {
        size_t old_prec = global_precision;
        global_precision = prec + 10;

        BigFloat sum("1");
        BigFloat term("1");
        size_t k = 1;

        while (true) {
            term = term / BigFloat(static_cast<int64_t>(k));
            if (term.is_zero()) break;
            sum = sum + term;
            k++;
        }

        global_precision = old_prec;
        sum._round(prec);
        return sum;
    }

private:
    void _parse_string(const std::string& s) {
        std::string str = s;
        // Trim whitespace
        str.erase(0, str.find_first_not_of(" \t\n\r"));
        str.erase(str.find_last_not_of(" \t\n\r") + 1);

        if (str.empty()) throw std::invalid_argument("Empty BigFloat string");

        if (str == "NaN" || str == "+NaN" || str == "-NaN") { is_nan = true; return; }
        if (str == "Infinity" || str == "+Infinity" || str == "inf") { is_inf = true; sign = 1; return; }
        if (str == "-Infinity" || str == "-inf") { is_inf = true; sign = -1; return; }

        if (str[0] == '+') { sign = 1; str = str.substr(1); }
        else if (str[0] == '-') { sign = -1; str = str.substr(1); }
        else { sign = 1; }

        size_t e_pos = str.find_first_of("eE");
        int64_t exp_offset = 0;
        if (e_pos != std::string::npos) {
            std::string exp_str = str.substr(e_pos + 1);
            str = str.substr(0, e_pos);
            exp_offset = std::stoll(exp_str);
        }

        size_t dot_pos = str.find('.');
        if (dot_pos != std::string::npos) {
            std::string int_part = str.substr(0, dot_pos);
            std::string frac_part = str.substr(dot_pos + 1);
            str = int_part + frac_part;
            exp_offset -= static_cast<int64_t>(frac_part.size());
        }

        // Lstrip zeros
        size_t first_nonzero = str.find_first_not_of('0');
        if (first_nonzero == std::string::npos) {
            sign = 0;
            mantissa = BigInt(0);
            exponent = 0;
            return;
        }

        str = str.substr(first_nonzero);
        mantissa = BigInt(str);
        exponent = exp_offset;
    }
};

inline std::ostream& operator<<(std::ostream& os, const BigFloat& bf) {
    os << bf.to_string();
    return os;
}

// Non-member scalar operators for interoperability
template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator+(const BigFloat& a, T b) { return a + BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator+(T a, const BigFloat& b) { return BigFloat(a) + b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator-(const BigFloat& a, T b) { return a - BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator-(T a, const BigFloat& b) { return BigFloat(a) - b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator*(const BigFloat& a, T b) { return a * BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator*(T a, const BigFloat& b) { return BigFloat(a) * b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator/(const BigFloat& a, T b) { return a / BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline BigFloat operator/(T a, const BigFloat& b) { return BigFloat(a) / b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator==(const BigFloat& a, T b) { return a == BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator==(T a, const BigFloat& b) { return BigFloat(a) == b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator!=(const BigFloat& a, T b) { return a != BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator!=(T a, const BigFloat& b) { return BigFloat(a) != b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator<(const BigFloat& a, T b) { return a < BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator<(T a, const BigFloat& b) { return BigFloat(a) < b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator<=(const BigFloat& a, T b) { return a <= BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator<=(T a, const BigFloat& b) { return BigFloat(a) <= b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator>(const BigFloat& a, T b) { return a > BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator>(T a, const BigFloat& b) { return BigFloat(a) > b; }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator>=(const BigFloat& a, T b) { return a >= BigFloat(b); }

template <typename T>
requires std::is_arithmetic_v<T>
inline bool operator>=(T a, const BigFloat& b) { return BigFloat(a) >= b; }

inline BigFloat abs(const BigFloat& in) {
    return in < 0 ? -in : in;
}

} // namespace math

namespace std {
inline math::BigFloat abs(const math::BigFloat& in) {
    return in < 0 ? -in : in;
}
} // namespace std

using math::abs;

namespace mat {
template <typename T> struct is_custom_scalar;
template <> struct is_custom_scalar<math::BigFloat> : std::true_type {};
}

#endif // BIGFLOAT_HPP
