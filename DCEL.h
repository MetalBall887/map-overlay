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

	halfEdge (Pt p, Pt q) : p (p), q (q) {}
};

struct Face {
	halfEdge* outer;
	std::vector <halfEdge*> inner;
	Edge* leftmost;
};

struct Vertex {
	double x, y;
	std::vector <halfEdge*> incident;

	Vertex (Pt p) : x (p.x), y (p.y) {}
};

struct DCEL {
	vector <Vertex> v;
	vector <Edge> e;
	vector <Face> f;
}