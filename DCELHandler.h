#include "DCEL.h"
using namespace std;

class DCELHandler {
	struct comp_triplets {
		bool operator () (const vector <double>& a, const vector <double>& b) const {
			for (int i = 0; i < 3; i++) {
				if (a[i] + EPS < b[i]) return true;
				if (abs (a[i] - b[i]) > EPS) return false;
			}
			return false;
		}
	};

	struct division_data {
		Pt x;
		int flag;
		halfEdge* ref;

		division_data (Pt x, int flag, halfEdge* ref) 
			: x (x), flag (flag), ref (ref) {}

		bool operator < (const division_data& b) {
			return (x == b.x ? flag < b.flag : x < b.x);
		}
	};

	static bool lex (const Edge& a, const Edge& b);
	static bool _radialComp (const Pt& a, const Pt& b);
	static inline bool up (const Pt& p);
	static bool radialComp (halfEdge* a, halfEdge* b);
	static bool intersect (const Edge& a, const Edge& b, double y, Pt& r);
	decltype (auto) next_it (set <Edge> :: iterator it, set <Edge>& s);
	decltype (auto) prev_it (set <Edge> :: iterator it, set <Edge>& s);
	vector <Edge> resolveOverlap (vector <Edge> e);
	pair < vector <Pt>, vector <Edge> > lineSegInt (vector <Edge> v);
	map <Pt, halfEdge*> findClosest (vector <halfEdge*> edges, vector <Pt> v);
	vector <double> get_line (Pt p, Pt q);
public:
	void report (halfEdge* a);
	void report (Edge a);
	DCEL construct (std::vector <Edge> e);
	void fill (DCEL& D, vector <Pt> p);
	DCEL merge (DCEL& A, DCEL& B);
};