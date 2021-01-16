
#ifndef D2_HH
#define D2_HH

#include <vector>

#include "ds++/Assert.hh"

class D2
{
public:
    D2() = default;
    explicit D2(unsigned N) : df(N), ddf(N, std::vector<double>(N)) {}
	D2(double d, unsigned N) : f(d), df(N, 0.0), ddf(N, std::vector<double>(N, 0.0)) {}

	D2(double d, unsigned i, unsigned N) : f(d), df(N, 0.0), ddf(N, std::vector<double>(N, 0.0)) 
    {
      Require(i<df.size());
      df[i] = 1.0;
    }

    unsigned size() const noexcept { return df.size(); }
    void swap(D2 &r);

    double y() const noexcept { return f; }
    double dydx(unsigned i) const { Require(i<df.size()); return df[i]; }
	double d2ydx2(unsigned i, unsigned j) const { Require(i<ddf.size()); Require(j<ddf[i].size()); return ddf[i][j]; } 

    friend D2 operator*(D2 const &, D2 const &);
    friend D2 operator*(D2 const &, double);
    friend D2 operator*(D2 &&, double);
    friend D2 operator*(double const a, D2 const &b){ return b*a; }
    D2 &operator*=(D2 const &);
    D2 &operator*=(double);

    friend D2 operator/(D2 const &, D2 const &);
    friend D2 operator/(D2 &&, D2 const &);
    friend D2 operator/(double, D2 const &);
    friend D2 operator/(double, D2 &&);
    D2 &operator/=(double const r){ return operator*=(1.0/r); }

    friend D2 operator+(D2 const &, D2 const &);
    friend D2 operator+(D2 const &, double);
    D2 &operator+=(D2 const &);

    friend D2 operator-(double, D2 const &);
    friend D2 operator-(D2 const &, D2 const &);
    friend D2 operator-(D2 &&, D2 const &);
    D2 &operator-=(D2 const &);

	friend D2 operator-(D2 const &);

    friend bool operator>(D2 const &a, double const b){ return a.f > b; }
    friend bool operator<(D2 const &a, D2 const &b){ return a.f < b.f; }
    friend bool operator<(D2 const &a, double const b){ return a.f < b; }
    friend bool operator<=(D2 const &a, D2 const &b){ return a.f <= b.f; }
    friend bool operator<=(D2 const &a, double const b){ return a.f <= b; }

    friend D2 log(D2 const &);

private:
  double f;
  std::vector<double> df;
  std::vector<std::vector<double>> ddf;
};

template<typename Real> unsigned size(Real const &);
template<> inline unsigned size(double const &){ return 1.0; }
template<> inline unsigned size(D2 const &x){ return x.size(); }

template<typename Real>
Real to_Real(double d, unsigned n);

template<>
inline double to_Real<double>(double d, unsigned){ return d; }

template<>
inline D2 to_Real<D2>(double d, unsigned n){ return D2(d, n); }

template<typename Real>
Real to_Real(double d, unsigned i, unsigned n);

template<>
inline double to_Real<double>(double d, unsigned, unsigned){ return d; }

template<>
inline D2 to_Real<D2>(double d, unsigned i, unsigned n){ return D2(d, i, n); }

#endif // D2_HH