#include <complex>
#include <iomanip>
#include <string>
#include <sstream>



template <typename T> class complexo:public std::complex<T>{ 
private:
   template<typename U> friend std::ostream &operator<<(std::ostream &stream,complexo<U> z);
   template<typename U> friend complexo operator *(double real,complexo<U> z){
        return complexo<U>(real*z.real(),real*z.imag());
   }
   template<typename U>friend complexo operator *(complexo<U> z,double real){
        return complexo<U>(real*z.real(),real*z.imag());
   }
   template<typename U> friend complexo operator +(double real,complexo<U> z){
          return complexo<U>(real+z.real(),z.imag());
   }       
   template<typename U> friend complexo operator +(complexo<U> z,double real){
          return complexo<U>(real+z.real(),z.imag());
   }
   template<typename U> friend complexo operator -(double real,complexo<U> z){
          return complexo<U>(real-z.real(),-z.imag());
   }       
   template<typename U> friend complexo operator -(complexo<U> z,double real){
          return complexo<U>(-real+z.real(),z.imag());
   }
   /*
   template<typename U> friend complexo operator /(double real,complexo<U> z){
             return complexo<U>(real/z.real(),real/z.imag());
   }
   template<typename U> friend complexo operator /(complexo<U> z,double real){
                 return complexo<U>(z.real()/real,z.imag()/real);
   }*/
public:
    complexo(T r=0,T i=0):std::complex<T>(r,i){}
    complexo(const std::complex<T> &c):std::complex<T>(c){}
    /*complexo<T> operator*(const complexo<T> &c){
        return complexo<T>(this->real()*c.real()-this->imag()*c.imag(),this->real()*c.imag()+this->imag()*c.real());
    }*//*
    complexo<T> operator+(const complexo<T> &c){
        return complexo<T>(this->real()+c.real(),this->imag()+c.imag());
    }*/
    complexo<T> operator-(const complexo<T> &c){
        return complexo<T>(this->real()-c.real(),this->imag()-c.imag());
    }
    complexo<T> operator*(const T s){
        return complexo<T>(this->real()*s,this->imag()*s);
    }
    complexo<T> operator/(const T s){
        return complexo<T>(this->real()/s,this->imag()/s);
    }
    /*
    complexo<T> operator/(const complexo<T> &c){
        return complexo<T>((this->real()*c.real()+this->imag()*c.imag())/(c.real()*c.real()+c.imag()*c.imag()),(this->imag()*c.real()-this->real()*c.imag())/(c.real()*c.real()+c.imag()*c.imag()));
    }*/
    /*complexo operator  *(complexo<T> z){
        return complexo<T>(this->real()*z.real()-this->imag()*z.imag(),this->real()*z.imag()+this->imag()*z.real());
    }*/
    /*
    complexo<T> operator^(const int n){
        complexo<T> c(1,0);
        for (int i = 0; i < n; ++i) {
            c = c*(*this);
        }
        return c;
    }*/
};

template <typename U>
std::string format_scalar(U val) {
    if (std::abs(val) < 1e-12) return "0";
    if (std::abs(val - std::round(val)) < 1e-9) {
        return std::to_string(static_cast<long long>(std::round(val)));
    }
    std::ostringstream ss;
    ss << std::defaultfloat << std::setprecision(5) << val;
    return ss.str();
}

template <typename U>
std::ostream &operator<<(std::ostream &stream, complexo<U> z) {
    U r = std::abs(z.real()) < 1e-12 ? U{0} : z.real();
    U im = std::abs(z.imag()) < 1e-12 ? U{0} : z.imag();

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
    //if (n < 0) return complexo<double>(0, std::sqrt(static_cast<double>(-n)));
    return complexo<double>(std::sqrt(static_cast<complexo<double>>(n)));
}

