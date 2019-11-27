#include "DCEL.h"

using namespace std;

Pt cur;
int ptr = 0;

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

struct Edge {
	Pt p, q;
	double A, B, C;
	halfEdge* origin;

	Edge () {
		p = Pt (0, 0);
		q = Pt (0, 0);
		if (tie (p.y, p.x) > tie (q.y, q.x)) swap (p, q);

		A = p.y - q.y;
		B = q.x - p.x;
		C = -p.x * A - q.x * B;
	}

	Edge (double x1, double y1, double x2, double y2) : origin (NULL) {
		p = Pt (x1, y1);
		q = Pt (x2, y2);
		if (tie (p.y, p.x) > tie (q.y, q.x)) swap (p, q);

		A = p.y - q.y;
		B = q.x - p.x;
		C = -p.x * A - q.x * B;
	}

	Edge (Pt p, Pt q, halfEdge* e) : p (p), q (q), origin (e) {
		if (tie (p.y, p.x) > tie (q.y, q.x)) swap (p, q);
		A = p.y - q.y;
		B = q.x - p.x;
		C = -p.x * A - q.x * B;
	}

	bool operator == (const Edge& b) const {
		return (p - b.p).dist () < EPS && (q - b.q).dist () < EPS;
	}

	double eval (double& x) {
		if (B) return -(A * x + C) / B;
		else return 1e9;
	}

	double rev (const double& y) const {
		if (A) return -(B * y + C) / A;
		else if (y == p.y) return cur.x;
		else if (p.y < y) return 1e9;
		else return -1e9;
	}

	bool operator < (const Edge& b) const {
		if (abs (rev (cur.y) - b.rev (cur.y)) > EPS) return rev (cur.y) < b.rev (cur.y);
		else {
			if (abs (rev (cur.y + 0.005) - b.rev (cur.y + 0.005)) < EPS) return p.y < b.p.y;
			if (cur.x <= rev (cur.y)) return rev (cur.y - 0.005) < b.rev (cur.y - 0.005);
			else return rev (cur.y + 0.005) < b.rev (cur.y + 0.005);
		}
	}
};

bool lex (const Edge& a, const Edge& b) {
	return tie (a.p, a.q) < tie (b.p, b.q);
}

bool _radialComp (const Pt& a, const Pt& b) {
	if (a.x < 0 && b.x >= 0) return false;
	if (a.x >= 0 && b.x < 0) return true;
	if (a.x == 0 && b.x == 0) return a.y < b.y;
	return a.cross (b) < 0;
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

pair < vector <Pt>, vector <Edge> > lineSegInt (vector <Edge> v) {
	vector <Pt> res;
	vector <Edge> res2;
	map < Pt, vector <Edge> > start, inter;
	map < Pt, set <Edge> > end;
	set <Pt> q;
	set <Edge> s;

	for (Edge a : v) {
		start[a.p].push_back (a);
		end[a.q].insert (a);
		q.insert (a.p.swap ());
		q.insert (a.q.swap ());
	}

	while (!q.empty ()) {
		auto x = *q.begin ();
		swap (x.x, x.y);

		cur = x;
		if (start.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			if (start[x].size () > 1) res.push_back (x);
			for (auto a : start[x]) {
				cout << '+' << ' ' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				auto it = s.insert (a).first;
				auto pr = prev_it (it, s), nx = next_it (it, s);
				Pt r;
				
				if (pr != s.end () && intersect (*it, *pr, x.y, r)) {
					//cout << "a!" << endl;
					inter[r].push_back (*it);
					inter[r].push_back (*pr);
					q.insert (r.swap ());
				}

				if (nx != s.end () && intersect (*nx, *it, x.y, r)) {
					//cout << "a!" << endl;
					inter[r].push_back (*it);
					inter[r].push_back (*nx);
					q.insert (r.swap ());
				}
			}
			start.erase (x);
		}

		if (inter.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			sort (inter[x].begin(), inter[x].end(), lex);
			int it = unique (inter[x].begin(), inter[x].end()) - inter[x].begin ();
			inter[x].resize (it);
			cout << it << endl;
			for (auto a : inter[x]) {
				cout << 'X' << ' ' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				if ((x - a.p).dist () >= EPS) res2.push_back (Edge (a.p, x, a.origin));
				s.erase (a);
				end[a.q].erase (a);
			}
			cur.x += 1e-3;

			for (auto a : inter[x]) {
				if ((a.q - x).dist () < EPS) continue;
				s.insert (Edge (x, a.q, a.origin));
				end[a.q].insert (Edge (x, a.q, a.origin));
			}

			inter.erase (x);
		}

		if (end.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			for (auto a : end[x]) {
				cout << '-' << ' ' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				if ((a.q - a.p).dist () < EPS) continue;
				res2.push_back (a);
				s.erase (a);
			}

			if (end[x].size () && s.size () > 1) {
				auto it = s.find (*end[x].begin ());
				auto pr = prev_it (it, s), nx = next_it (it, s);
				Pt r;

				if (nx != s.end () && pr != s.end () && intersect (*nx, *pr, x.y + EPS, r)) {
					inter[r].push_back (*pr);
					inter[r].push_back (*nx);
					q.insert (r.swap ());
				}	
			}
			end.erase (x);
		}

		for (auto b : s)
			cout << "Set: " << b.p.x << ' ' << b.p.y << ' ' << b.q.x << ' ' << b.q.y << endl;
		cout << endl;

		q.erase (q.begin ());
	}

	sort (res.begin(), res.end());
	int it = unique (res.begin(), res.end()) - res.begin ();
	res.resize (it);

	for (auto a : res2)
		cout << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;

	return make_pair (res, res2);
}

map <Pt, halfEdge*> findClosest (vector <Edge> e, vector <Pt> v) {
	set <Pt> q;
	set <Edge> s;
	map < Pt, vector <Edge> > start, end;
	map < Pt, vector <Pt> > points;
	map <Pt, halfEdge*> res;
	Edge border (-1e5, -1e9, -1e5, 1e9);
	s.insert (border);

	for (auto a : e) {
		start[a.p].push_back (a);
		end[a.q].push_back (a);
		q.insert (a.p.swap ());
		q.insert (a.q.swap ());
	}

	for (auto a : v) {
		points[a].push_back (a);
		q.insert (a.swap ());
	}

	for (Pt x : q) {
		x = x.swap ();
		cur = x;
		cout << "Sweeping line at: " << x.x << ' ' << x.y << endl;
		if (start.count (x)) {
			for (auto a : start[x]) {
				cout << '+' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				s.insert (a);
			}
		}

		if (points.count (x)) {
			for (auto a : points[x]) {
				auto it = s.lower_bound (Edge (a.x - EPS, -1e9, a.x - EPS, 1e9));
				it--;
				res[a] = it -> origin;
			}
		}

		if (end.count (x)) {
			for (auto a : end[x]) {
				cout << '-' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				s.erase (a);
			}
		}

		for (auto b : s)
			cout << "Set: " << b.p.x << ' ' << b.p.y << ' ' << b.q.x << ' ' << b.q.y << endl;
		cout << endl;
	}

	return res;
}