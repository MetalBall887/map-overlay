#include "DCELHandler.h"

struct Face;
struct halfEdge;
struct Vertex;


Edge::Edge ();
Edge::Edge (double x1, double y1, double x2, double y2);
Edge::Edge (Pt p, Pt q, halfEdge* e);
Edge::Edge (const halfEdge& e);

double Edge::rev (const double& y) const;

bool Edge:: operator== (const Edge& b) const;
bool Edge:: operator< (const Edge& b) const;

Edge::halfEdge (Pt p, Pt q) 
	: p (p), q (q), incidentFace (NULL), to_a (NULL), to_b (NULL) {}
	
Edge::halfEdge () : incidentFace (NULL) {}

Face::~Face () {
	inner.clear ();
}

Vertex::Vertex (Pt p) : p (p) {}

bool Vertex::operator < (const Vertex& b) const {
	return p < b.p;
}

Vertex::~Vertex () {
	incident.clear ();
}

struct DCEL {
	std::vector <Vertex*> v;
	std::vector <halfEdge*> e;
	std::vector <Face*> f;
	Face* outer_bound;

	~DCEL () {
		for (auto a : v) delete a;
		for (auto a : e) delete a;
		for (auto a : f) delete a;
		v.clear (), e.clear (), f.clear ();
	}
};