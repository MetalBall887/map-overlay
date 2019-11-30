#include "Point.h"

typedef Point <double> Pt;

struct Face;
struct halfEdge;
struct Vertex;

struct halfEdge {
	Pt p, q;
	halfEdge *twin, *next, *prev;
	Face* incidentFace;
	Vertex* origin;

	halfEdge (Pt p, Pt q) : p (p), q (q), incidentFace (NULL) {}

	halfEdge () : incidentFace (NULL) {}
};

struct Face {
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

	~DCEL () {
		for (auto a : v) delete a;
		for (auto a : e) delete a;
		for (auto a : f) delete a;
		v.clear (), e.clear (), f.clear ();
	}
};