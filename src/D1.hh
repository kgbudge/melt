
#ifndef D1_HH
#define D1_HH

#include <vector>

#include "ds++/Assert.hh"

class D1
{
public:
    D1() = default;
    explicit D1(unsigned N) : df(N) {}
	D1(double d, unsigned N) : f(d), df(N, 0.0) {}

	D1(double d, unsigned i, unsigned N) : f(d), df(N, 0.0)
    {
      Require(i<df.size());
      df[i] = 1.0;
    }

    unsigned size() const noexcept { return df.size(); }
    void swap(D1 &r);

    void dydx(unsigned i, double x){ Require(i<df.size()); df[i] = x; }

    explicit operator double() const { return f; }

    double y() const noexcept { return f; }
    double dydx(unsigned i) const { Require(i<df.size()); return df[i]; }

    friend D1 operator*(D1 const &, D1 const &);
    friend D1 operator*(D1 const &, double);
    friend D1 operator*(D1 &&, double);
    friend D1 operator*(double const a, D1 const &b){ return b*a; }
    D1 &operator*=(D1 const &);
    D1 &operator*=(double);

    friend D1 operator/(D1 const &, D1 const &);
    friend D1 operator/(D1 &&, D1 const &);
    friend D1 operator/(double, D1 const &);
    friend D1 operator/(double, D1 &&);
    D1 &operator/=(D1 const &);
    D1 &operator/=(double const r){ return operator*=(1.0/r); }

    friend D1 operator+(D1 const &, D1 const &);
    friend D1 operator+(D1 const &, double);
    D1 &operator+=(double);
    D1 &operator+=(D1 const &);

    friend D1 operator-(double, D1 const &);
    friend D1 operator-(D1 const &, D1 const &);
    friend D1 operator-(D1 &&, D1 const &);
    D1 &operator-=(double);
    D1 &operator-=(D1 const &);

	friend D1 operator-(D1 const &);

    friend bool operator>(D1 const &a, double const b){ return a.f > b; }
    friend bool operator<(D1 const &a, D1 const &b){ return a.f < b.f; }
    friend bool operator<(D1 const &a, double const b){ return a.f < b; }
    friend bool operator<=(D1 const &a, D1 const &b){ return a.f <= b.f; }
    friend bool operator<=(D1 const &a, double const b){ return a.f <= b; }

    friend D1 log(D1 const &);

private:
  double f;
  std::vector<double> df;
};

template<typename Real> unsigned size(Real const &);
template<> inline unsigned size(double const &){ return 1.0; }
template<> inline unsigned size(D1 const &x){ return x.size(); }

template<typename Real>
Real to_Real(double d, unsigned n);

template<>
inline double to_Real<double>(double d, unsigned){ return d; }

template<>
inline D1 to_Real<D1>(double d, unsigned n){ return D1(d, n); }

template<typename Real>
Real to_Real(double d, unsigned i, unsigned n);

template<>
inline double to_Real<double>(double d, unsigned, unsigned){ return d; }

template<>
inline D1 to_Real<D1>(double d, unsigned i, unsigned n){ return D1(d, i, n); }



template<typename Real>
double value(Real const &);

template<>
inline double value(double const &x) { return x; }

template<>
inline double value(D1 const &x) { return x.y(); }



template<typename Real>
double dydx(Real const &, int);

template<>
inline double dydx(D1 const &x, int i){ return x.dydx(i);}

template<>
inline double dydx(double const &x, int ) { return 0.0; }

#endif // D1_HH