#include "Point.h"

struct Face;
struct halfEdge;
struct Vertex;

struct Edge {
	Pt p, q;
	double A, B, C;
	halfEdge* origin;

	Edge ();
	Edge (double x1, double y1, double x2, double y2);
	Edge (Pt p, Pt q, halfEdge* e);
	Edge (const halfEdge& e);

	double rev (const double& y) const;

	bool operator == (const Edge& b) const;
	bool operator < (const Edge& b) const;
};

struct halfEdge {
	Pt p, q;
	halfEdge *twin, *next, *prev;
	Face* incidentFace;
	Vertex* origin;
	Face *to_a, *to_b;

	halfEdge (Pt p, Pt q) 
		: p (p), q (q), incidentFace (NULL), to_a (NULL), to_b (NULL) {}
	
	halfEdge () : incidentFace (NULL) {}
};

struct Face {
	bool painted;
	Face *to_a, *to_b;
	halfEdge* outer;
	std::vector <halfEdge*> inner;
	halfEdge* leftmost;

	~Face () {
		inner.clear ();
	}
};

struct Vertex {
	Pt p;
	std::vector <halfEdge*> incident;
	Vertex (Pt p) : p (p) {}

	bool operator < (const Vertex& b) const {
		return p < b.p;
	}

	~Vertex () {
		incident.clear ();
	}
};

struct DCEL {
	std::vector <Vertex*> v;
	std::vector <halfEdge*> e;
	std::vector <Face*> f;
	Face* outer_bound;

	DCEL () {}

	~DCEL () {
		for (auto a : v) delete a;
		for (auto a : e) delete a;
		for (auto a : f) delete a;
		v.clear (), e.clear (), f.clear ();
	}

	DCEL (const DCEL&& other) {
		v = other.v;
		e = other.e;
		f = other.f;
		outer_bound = other.outer_bound;
	}

	DCEL (const DCEL& other);
	DCEL& operator= (const DCEL& other) {
		if (this == &other) return *this;
		DCEL tmp (other);
		std::swap (v, tmp.v);
		std::swap (e, tmp.e);
		std::swap (f, tmp.f);
		outer_bound = other.outer_bound;
		return *this;
	}
};