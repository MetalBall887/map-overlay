#include "DCEL.h"

using namespace std;

Pt cur;

void report (halfEdge* a) {
	cout << a -> p.x << ' ' << a -> p.y << ' ' << a -> q.x << ' ' << a -> q.y << endl;
}

void report (Edge a) {
	cout << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
}

int segmentIntersection (const Pt& s1, const Pt& e1,
		const Pt& s2, const Pt& e2, Pt& r1, Pt& r2) {
	if (e1 == s1) {
		if (e2 == s2) {
			if (e1 == e2) { r1 = e1; return 1; } //all equal
			else return 0; //different point segments
		} else return segmentIntersection(s2,e2,s1,e1,r1,r2);//swap
	}
	//segment directions and separation
	Pt v1 = e1 - s1, v2 = e2 - s2, d = s2 - s1;
	auto a = v1.cross (v2), a1 = v1.cross (d), a2 = v2.cross (d);
	if (a == 0) { //if parallel
		auto b1 = s1.dot (v1), c1 = e1.dot (v1),
		     b2 = s2.dot (v1), c2 = e2.dot (v1);
		if (a1 || a2 || max(b1, min (b2,c2)) > min (c1, max (b2,c2)))
			return 0;
		r1 = min (b2,c2) < b1 ? s1 : (b2 < c2 ? s2 : e2);
		r2 = max (b2, c2) > c1 ? e1 : (b2 > c2 ? s2 : e2);
		return 2 - (r1 == r2);
	}
	if (a < 0) { a = -a; a1 = -a1; a2 = -a2; }
	if (0 < a1 || a < -a1 || 0 < a2 || a < -a2)
		return 0;
	r1 = s1 - v1 * a2 / a;
	return 1;
}

Edge::Edge () {
	p = Pt (0, 0);
	q = Pt (0, 0);
	if (p > q) swap (p, q);

	A = p.y - q.y;
	B = q.x - p.x;
	C = -p.x * A - p.y * B;
}

Edge::Edge (double x1, double y1, double x2, double y2) : origin (NULL) {
	p = Pt (x1, y1);
	q = Pt (x2, y2);
	if (p > q) swap (p, q);

	A = p.y - q.y;
	B = q.x - p.x;
	C = -p.x * A - p.y * B;
}

Edge::Edge (Pt a, Pt b, halfEdge* e) : p (a), q (b), origin (e) {
	if (p > q) swap (p, q);
	A = p.y - q.y;
	B = q.x - p.x;
	C = -p.x * A - p.y * B;
}

Edge::Edge (const halfEdge& e) : p (e.p), q (e.q), origin (NULL) {
	if (p > q) swap (p, q);

	A = p.y - q.y;
	B = q.x - p.x;
	C = -p.x * A - p.y * B;	
}

bool Edge::operator == (const Edge& b) const {
	return (p - b.p).dist () < EPS && (q - b.q).dist () < EPS;
}

double Edge::rev (const double& y) const {
	if (A) return -(B * y + C) / A;
	else if (y == p.y) return p.x;
	else if (y < p.y) return -1e9;
	else return 1e9;
}

bool Edge::operator < (const Edge& b) const {
	Pt cr, cr2;
	segmentIntersection (p, q, b.p, b.q, cr, cr2);
	if (abs (rev (cur.y) - b.rev (cur.y)) > EPS) return rev (cur.y) < b.rev (cur.y);
	else {
		if (abs (rev (cur.y + 0.005) - b.rev (cur.y + 0.005)) < EPS) return p < b.p;
		if ((cr - p).dist () < EPS && (cr - b.p).dist () < EPS) return rev (cur.y + 0.005) < b.rev (cur.y + 0.005);
		if ((cr - p).dist () < EPS) return true;
		if ((cr - b.p).dist () < EPS) return false;
		return rev (cur.y - 0.005) < b.rev (cur.y - 0.005);
	}
}

bool lex (const Edge& a, const Edge& b) {
	return tie (a.p, a.q) < tie (b.p, b.q);
}

inline bool up (const Pt& p) {
  return p.y > 0 or (p.y == 0 and p.x >= 0);
}

bool _radialComp (const Pt& a, const Pt& b) {
	return up(a) == up(b) ? a.x * b.y < a.y * b.x : up(a) > up(b);
}

bool intersect (const Edge& a, const Edge& b, double y, Pt& r) {
	Pt r1, r2;
	int res = segmentIntersection (a.p, a.q, b.p, b.q, r1, r2);
	r = r1;
	if (res == 1 && r.y >= y) return true;
	return false;
}

decltype (auto) next_it (set <Edge> :: iterator it, set <Edge>& s) {
	it++;
	return it;
}

decltype (auto) prev_it (set <Edge> :: iterator it, set <Edge>& s) {
	if (it == s.begin ()) return s.end ();
	it--;
	return it;
}

struct comp_triplets {
	bool operator () (const vector <double>& a, const vector <double>& b) const {
		for (int i = 0; i < 3; i++) {
			if (a[i] + EPS < b[i]) return true;
			if (abs (a[i] - b[i]) > EPS) return false;
		}
		return false;
	}
};

vector <Edge> resolveOverlap (vector <Edge> e) {
	vector <Edge> ans;
	map <vector <double>, vector<pair <Pt, bool>>, comp_triplets> m;
	int cnt = 0;

	for (auto x : e) {
		double A = x.A, B = x.B, C = x.C;
		if (A) B /= A, C /= A, A /= A;
		else C /= B, B /= B;
		vector <double> v {A, B, C};
		m[v].push_back ({x.p, true});
		m[v].push_back ({x.q, false});
	}

	for (auto a : m) {
		auto& v = a.second;
		sort (v.begin(), v.end());
		for (int i = 0; i < v.size (); i++) {
			if (v[i].second) cnt++;
			else cnt--;
			if (cnt) {
				if (v[i].first != v[i+1].first) {
					ans.push_back (Edge (v[i].first, v[i+1].first, NULL));
					report (ans.back ());
				}
			}
		}
	}

	return ans;
}

pair < vector <Pt>, vector <Edge> > lineSegInt (vector <Edge> v) {
	v = resolveOverlap (v);
	vector <Pt> res;
	vector <Edge> res2;
	map < Pt, vector <Edge> > start, inter, end;
	set <Pt> q;
	set <Edge> s;

	for (Edge a : v) {
		start[a.p].push_back (a);
		end[a.q].push_back (a);
		q.insert (a.p);
		q.insert (a.q);
	}

	for (auto& a : start) {
		sort (a.second.begin(), a.second.end(), lex);
	}

	for (auto& a : end) {
		sort (a.second.begin(), a.second.end(), lex);
	}

	while (!q.empty ()) {
		auto x = *q.begin ();

		cur = x;
		if (start.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			if (start[x].size () > 1) res.push_back (x);
			for (auto a : start[x]) {
				auto it = s.insert (a).first;
				auto pr = prev_it (it, s), nx = next_it (it, s);
				Pt r;

				if (pr != s.end () && intersect (*it, *pr, x.y, r)) {
					inter[r].push_back (*it);
					inter[r].push_back (*pr);
					q.insert (r);
				}

				if (nx != s.end () && intersect (*nx, *it, x.y, r)) {
					inter[r].push_back (*it);
					inter[r].push_back (*nx);
					q.insert (r);
				}
			}
			start.erase (x);
		}

		if (inter.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			sort (inter[x].begin(), inter[x].end(), lex);
			int sz = unique (inter[x].begin(), inter[x].end()) - inter[x].begin ();
			inter[x].resize (sz);
			for (auto a : inter[x]) {
				if (s.count (a) && (x - a.p).dist () >= EPS) res2.push_back (Edge (a.p, x, a.origin));
				s.erase (a);
			}

			for (auto a : inter[x]) {
				if ((a.q - x).dist () < EPS) continue;
				auto it = s.insert (Edge (x, a.q, a.origin)).first;

				auto pr = prev_it (it, s), nx = next_it (it, s);
				Pt r;

				if (pr != s.end () && intersect (*it, *pr, x.y, r)) {
					inter[r].push_back (*it);
					inter[r].push_back (*pr);
					q.insert (r);
				}

				if (nx != s.end () && intersect (*nx, *it, x.y, r)) {
					inter[r].push_back (*it);
					inter[r].push_back (*nx);
					q.insert (r);
				}

				end[a.q].push_back (Edge (x, a.q, a.origin));
			}

			inter.erase (x);
		}

		if (end.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			Edge ex;
			for (auto a : end[x]) {
				ex = a;
				if (s.find (a) != s.end () && (a.q - a.p).dist () >= EPS) res2.push_back (a);
				s.erase (a);
			}
			auto it = s.lower_bound (ex);
			if (it != s.end () && end[x].size () && s.size () > 1) {
				auto pr = prev_it (it, s);
				Pt r;

				if (pr != s.end () && intersect (*pr, *it, x.y + EPS, r)) {
					inter[r].push_back (*pr);
					inter[r].push_back (*it);
					q.insert (r);
				}
			}
			end.erase (x);
		}
		q.erase (q.begin ());
	}
	sort (res.begin(), res.end());
	int it = unique (res.begin(), res.end()) - res.begin ();
	res.resize (it);

	return make_pair (res, res2);
}

map <Pt, halfEdge*> findClosest (vector <halfEdge*> edges, vector <Pt> v) {
	vector <Edge> e;
	set <double> q;
	set <Edge> s;
	map < double, vector <Edge> > start, end;
	map < double, vector <Pt> > points;
	map <Pt, halfEdge*> res;
	Edge border (-1e6, -1e9, -1e6, 1e9);
	s.insert (border);

	for (auto a : edges) {
		e.push_back (Edge (a -> p, a -> q, a));
		assert (e.back ().p < e.back ().q);
	}

	for (auto a : e) {
		start[a.p.y].push_back (a);
		end[a.q.y].push_back (a);
		q.insert (a.p.y);
		q.insert (a.q.y);
	}

	for (auto a : v) {
		points[a.y].push_back (a);
		q.insert (a.y);
	}

	for (auto x : q) {
		cur = Pt (0, x);
		if (start.count (x)) {
			for (auto a : start[x]) {
				s.insert (a);
			}
		}

		if (points.count (x)) {
			for (auto a : points[x]) {
				auto it = s.lower_bound (Edge (a.x, a.y, a.x - 1, a.y + EPS));
				it--; 
				res[a] = it -> origin;
			}
		}

		if (end.count (x)) {
			for (auto a : end[x]) {
				s.erase (a);
			}
		}
	}

	q.clear (), s.clear (), start.clear (), end.clear (), points.clear ();

	return res;
}
