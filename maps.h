#include <iostream>
#include "line_segment.h"

using namespace std;

bool radialComp (halfEdge* a, halfEdge* b) {
	return _radialComp (a -> p - a -> q, b -> p - b -> q);
}

void report (halfEdge* a) {
	cout << a -> p.x << ' ' << a -> p.y << ' ' << a -> q.x << ' ' << a -> q.y << endl;
}

DCEL construct (std::vector <Edge> e) {
	auto res = lineSegInt (e);

	DCEL D;

	map <Pt, std::vector <halfEdge*> > origins;

	for (auto a : res.second) {
		halfEdge *fst = new halfEdge (a.p, a.q), *snd  = new halfEdge (a.q, a.p);
		fst -> twin = snd;
		snd -> twin = fst;

		origins[a.q].push_back (fst);
		origins[a.p].push_back (snd);

		D.e.push_back (fst);
		D.e.push_back (snd);
		cout << "-> ";
		report (fst);
		cout << "<- ";
		report (snd);
	}

	for (auto a : res.first) {
		Vertex* vk = new Vertex (a);
		D.v.push_back (vk);
		for (auto b : origins[a]) {
			vk -> incident.push_back (b);
			b -> origin = vk;
		}

		auto& hedgehog = vk -> incident;
		int n = hedgehog.size ();

		sort (hedgehog.begin (), hedgehog.end (), radialComp);

		for (int i = 0; i < n - 1; i++) {
			hedgehog[i] -> next = hedgehog[i + 1] -> twin;
			hedgehog[i + 1] -> twin -> prev = hedgehog[i];
		}
		hedgehog[n - 1] -> next = hedgehog[0] -> twin ;
		hedgehog[0] -> twin -> prev = hedgehog[n - 1];
	}

	set <halfEdge*> went;

	vector <Edge> edgeClosest;
	vector <Pt> tryClosest;
	map <Pt, Face*> waitClosest;

	for (halfEdge* a : D.e) {
		if (!a -> incidentFace) {
			Face* f = new Face ();
			f -> painted = false;
			a -> incidentFace = f;
			D.f.push_back (f);
			f -> leftmost = a;

			halfEdge* x = a -> next;

			while (x != a) {
				x -> incidentFace = a -> incidentFace;
				x = x -> next;
				if (f -> leftmost -> origin -> p > x -> origin -> p) {
					f -> leftmost = x;
				}
			}

			Pt pa = f -> leftmost -> twin -> origin -> p;
			Pt pb = f -> leftmost -> origin -> p;
			Pt pc = f -> leftmost -> next -> origin -> p;


			if ((pa - pb).cross (pc - pb) >= 0 || f -> leftmost -> next == f -> leftmost -> twin) {
				if (a -> p < a -> q) edgeClosest.push_back (Edge (a -> p, a -> q, a));
				x = a -> next;
				report (a);
				while (x != a) {
					report (x);
					if (x -> p > x -> q) edgeClosest.push_back (Edge (x -> q, x -> p, a));
					x = x -> next;
				}

				tryClosest.push_back (f -> leftmost -> origin -> p);
				waitClosest[f -> leftmost -> origin -> p] = f;
			}
			else {
				if (a -> p > a -> q) edgeClosest.push_back (Edge (a -> q, a -> p, a));
				x = a -> next;
				while (x != a) {
					if (x -> p > x -> q) edgeClosest.push_back (Edge (x -> q, x -> p, a));
					x = x -> next;
				}
			}
		}
	}

	auto mp = findClosest (edgeClosest, tryClosest);
	map <Face*, vector <Face*> > g;
	set <Face*> s;

	for (auto it : mp) {
		if (it.second) {
			g[it.second -> incidentFace].push_back (waitClosest[it.first]);
			s.insert (waitClosest[it.first]);
			cout << it.second -> incidentFace << ' ' << waitClosest[it.first] << endl;
		}
	}

	for (auto a : mp) {
		cout << a.first.x << ' ' << a.first.y << '+';
		if (a.second) report (a.second);
		else cout << "Outer bound";
	}

	queue < pair <Face*, Face*> > q;
	vector <Face*> new_f;
	for (Face* f : D.f) {
		if (!s.count (f)) {
			q.push ({f, f});
			new_f.push_back (f);
		}
	}

	while (!q.empty ()) {
		Face* x = q.front ().first;
		Face* p = q.front ().second;
		q.pop ();

		if (x != p) p -> inner.push_back (x -> leftmost);
		else { x -> outer = x -> leftmost; }

		auto a = x -> leftmost;
		a -> incidentFace = p;
		auto b = a -> next;
		while (b != a) { b -> incidentFace = p; b = b -> next;}

		for (auto to : g[x]) {
			q.push ({to, p});
		}
	}

	D.f.clear ();

	D.f = new_f;



	return D;
}

void fill (DCEL& D, vector <Pt> p) {
	vector <Edge> edges;
	for (Face* f : D.f) {
		auto a = f -> outer;

		if (a -> q < a -> p) edges.push_back (Edge (a -> q, a -> p, a));
		auto x = a -> next;

		while (x != a) {
			if (x -> q > x -> p) edges.push_back (Edge (x -> q, x -> p, x));
			x = x -> next;
		}

		for (auto a : f -> inner) {
			if (a -> p < a -> q) edges.push_back (Edge (a -> p, a -> q, a));
			auto x = a -> next;

			while (x != a) {
				if (x -> p > x -> q) edges.push_back (Edge (x -> p, x -> q, x));
				x = x -> next;
			}
		}
	}

	auto res = findClosest (edges, p);

	for (auto a : res) {
		a.second -> incidentFace -> painted = true;
	}

	for (auto f : D.f) {
		if (f -> painted == false) continue;

		auto a = f -> outer;

		report (a);
		auto x = a -> next;

		while (x != a) {
			report (x);
			x = x -> next;
		}

		cout << endl;
	}
}