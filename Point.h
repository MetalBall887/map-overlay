#include <iostream>
#include <cstdlib>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <set>
#include <cmath>
#include <map>

#define EPS 1e-6


template <class T>
struct Point {
	typedef Point P;
	T x, y;
	explicit Point (T x = 0, T y = 0) : x (x), y (y) {}
	bool operator < (P p) const { return ((y == p.y ? x < p.x : y < p.y) && p != *this); }
	bool operator > (P p) const { return p < *this; }
	bool operator == (P p) const { return (*this - p).dist () < EPS; }
	bool operator != (P p) const { return (*this - p).dist () >= EPS; }
	P operator + (P p) const { return P (x + p.x, y + p.y); }
	P operator - (P p) const { return P (x - p.x, y - p.y); }
	P operator * (T d) const { return P (x * d, y * d); }
	P operator / (T d) const { return P (x / d, y / d); }
	T dot (P p) const { return x * p.x + y * p.y; }
	T cross (P p) const { return x * p.y - y * p.x; }
	T cross (P a, P b) const { return (a - *this).cross (b - *this); }
	T dist2 () const { return x * x + y * y; }
	double dist () const { return sqrt ((double) dist2 ()); }
	double angle () const { return atan2 (y, x); }
	P unit () const { return *this / dist (); }
	P perp () const { return P (-y, x); }
	P swap () const { return P (y, x); }
	P normal () const { return perp ().unit (); }
	P rotate (double a) const {
		return P (x * cos (a) - y * sin (a), x * sin (a) + y * cos (a)); }
};
