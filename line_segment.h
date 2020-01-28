#include "DCEL.h"

using namespace std;

Pt cur;
bool ins;
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
		C = -p.x * A - p.y * B;
	}

	Edge (double x1, double y1, double x2, double y2) : origin (NULL) {
		p = Pt (x1, y1);
		q = Pt (x2, y2);
		if (tie (p.y, p.x) > tie (q.y, q.x)) swap (p, q);

		A = p.y - q.y;
		B = q.x - p.x;
		C = -p.x * A - p.y * B;
	}

	Edge (Pt p, Pt q, halfEdge* e) : p (p), q (q), origin (e) {
		if (tie (p.y, p.x) > tie (q.y, q.x)) swap (p, q);
		A = p.y - q.y;
		B = q.x - p.x;
		C = -p.x * A - p.y * B;
	}

	Edge (halfEdge& e) : p (e.p), q (e.q), origin (&e) {
		if (tie (p.y, p.x) > tie (q.y, q.x)) {
			origin = e.twin;
			swap (q, p);
		}
	}

	bool operator == (const Edge& b) const {
		return (p - b.p).dist () < EPS && (q - b.q).dist () < EPS;
	}

	double rev (const double& y) const {
		if (A) return -(B * y + C) / A;
		else if (y == p.y) return p.x;
		else if (y < p.y) return -1e9;
		else return 1e9;
	}

	bool operator < (const Edge& b) const {
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

vector <Edge> resolveOverlap (vector <Edge> e) {
	vector <Edge> ans;
	map <vector <double>, vector<pair <Pt, bool>>> m;
	int cnt = 0;

	for (auto x : e) {
		double A = x.A, B = x.B, C = x.C;
		if (A) A /= A, B /= A, C /= A;
		else B /= B, C /= B;
		vector <double> v {A, B, C};
		cout << v[0] << ' ' << v[1] << ' ' << v[2] << endl;
		m[v].push_back ({x.p, true});
		m[v].push_back ({x.q, false});
	}

	for (auto a : m) {
		cout << a.first[0] << ' ' << a.first[1] << ' ' << a.first[2] << endl;
		auto& v = a.second;
		sort (v.begin(), v.end());
		for (int i = 0; i < v.size (); i++) {
			if (v[i].second) cnt++;
			else cnt--;
			if (cnt) {
				if (v[i].first != v[i+1].first)
					ans.push_back (Edge (v[i].first, v[i+1].first, NULL));
				cout << v[i].first.x << ' ' << v[i].first.y << ' ' << v[i+1].first.x << ' ' << v[i+1].first.y << endl;
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

	while (!q.empty ()) {
		auto x = *q.begin ();

		cur = x;
		if (start.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			if (start[x].size () > 1) res.push_back (x);
			for (auto a : start[x]) {
				cout << '+' << ' ' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << ' ' << cur.x << ' ' << cur.y << endl;
				auto it = s.insert (a).first;
				auto pr = prev_it (it, s), nx = next_it (it, s);
				Pt r;

				if (pr != s.end () && intersect (*it, *pr, x.y, r)) {
					//cout << "a!" << endl;
					inter[r].push_back (*it);
					inter[r].push_back (*pr);
					q.insert (r);
				}

				if (nx != s.end () && intersect (*nx, *it, x.y, r)) {
					//cout << "a!" << endl;
					inter[r].push_back (*it);
					inter[r].push_back (*nx);
					q.insert (r);
				}
			}
			start.erase (x);
		}

		for (auto b : s)
			cout << "Set: " << b.p.x << ' ' << b.p.y << ' ' << b.q.x << ' ' << b.q.y << endl;
		cout << endl;

		if (inter.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			sort (inter[x].begin(), inter[x].end(), lex);
			int sz = unique (inter[x].begin(), inter[x].end()) - inter[x].begin ();
			inter[x].resize (sz);
			for (auto a : inter[x]) {
				if (s.find (a) == s.end ()) cout << "NOT IN SET\n";
				cout << 'X' << ' ' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				if (s.count (a) && (x - a.p).dist () >= EPS) res2.push_back (Edge (a.p, x, a.origin));
				s.erase (a);
				cout << s.size () << endl;
			}

			for (auto a : inter[x]) {
				if ((a.q - x).dist () < EPS) continue;
				auto it = s.insert (Edge (x, a.q, a.origin)).first;

				auto pr = prev_it (it, s), nx = next_it (it, s);
				Pt r;

				if (pr != s.end () && intersect (*it, *pr, x.y, r)) {
					//cout << "a!" << endl;
					inter[r].push_back (*it);
					inter[r].push_back (*pr);
					q.insert (r);
				}

				if (nx != s.end () && intersect (*nx, *it, x.y, r)) {
					//cout << "a!" << endl;
					inter[r].push_back (*it);
					inter[r].push_back (*nx);
					q.insert (r);
				}

				end[a.q].push_back (Edge (x, a.q, a.origin));
			}
			ins = false;

			inter.erase (x);
		}

		for (auto b : s)
			cout << "Set: " << b.p.x << ' ' << b.p.y << ' ' << b.q.x << ' ' << b.q.y << endl;
		cout << endl;

		if (end.count (x)) {
			if (res.empty () || res.back () != x) res.push_back (x);
			Edge ex;
			for (auto a : end[x]) {
				ex = a;
				cout << '-' << ' ' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
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
		cout << cur.x << ' ' << cur.y << ' ' << "AE \n";
		cout << q.size () << endl;
		for (auto b : s)
			cout << "Set: " << b.p.x << ' ' << b.p.y << ' ' << b.q.x << ' ' << b.q.y << endl;
		q.erase (q.begin ());
	}
	sort (res.begin(), res.end());
	int it = unique (res.begin(), res.end()) - res.begin ();
	res.resize (it);

	return make_pair (res, res2);
}

map <Pt, halfEdge*> findClosest (vector <Edge> e, vector <Pt> v) {
	set <double> q;
	set <Edge> s;
	map < double, vector <Edge> > start, end;
	map < double, vector <Pt> > points;
	map <Pt, halfEdge*> res;
	Edge border (-1e6, -1e9, -1e6, 1e9);
	s.insert (border);

	for (auto a : e) {
		start[a.p.y].push_back (a);
		end[a.q.y].push_back (a);
		q.insert (a.p.y);
		q.insert (a.q.y);
		cout << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
	}
	cout << endl;

	for (auto a : v) {
		cout << a.x << ' ' << a.y << endl;
		points[a.y].push_back (a);
		q.insert (a.y);
	}

	for (auto x : q) {
		cur = Pt (0, x);
		cout << "Sweeping line at: " << x << endl;
		if (start.count (x)) {
			for (auto a : start[x]) {
				cout << '+' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				s.insert (a);
			}
		}

		if (points.count (x)) {
			for (auto a : points[x]) {
				auto it = s.lower_bound (Edge (a.x, a.y, a.x - 1, a.y + EPS));
				it--;
				cout << '/' << a.x << ' ' << a.y << ' ' << it -> p.x << ' ' << it -> p.y << ' ' << it -> q.x << ' ' << it -> q.y << endl;
				res[a] = it -> origin;
			}
		}

		for (auto b : s)
			cout << "Set: " << b.p.x << ' ' << b.p.y << ' ' << b.q.x << ' ' << b.q.y << endl;
		cout << endl;

		if (end.count (x)) {
			for (auto a : end[x]) {
				cout << '-' << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl;
				s.erase (a);
			}
		}
	}

	q.clear (), s.clear (), start.clear (), end.clear (), points.clear ();

	return res;
}
