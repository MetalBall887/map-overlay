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

		fst -> origin = new Vertex (a.q);
		snd -> origin = new Vertex (a.p);

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
		Vertex vk (a);
		for (auto b : origins[a]) {
			vk.incident.push_back (b);
		}

		auto& hedgehog = vk.incident;
		int n = hedgehog.size ();

		sort (hedgehog.begin (), hedgehog.end (), radialComp);

		for (int i = 0; i < n - 1; i++) {
			hedgehog[i] -> next = hedgehog[i + 1] -> twin;
			hedgehog[i + 1] -> twin -> prev = hedgehog[i];

			report (hedgehog[i]);
			report (hedgehog[i + 1] -> twin);
		}
		hedgehog[n - 1] -> next = hedgehog[0] -> twin ;
		hedgehog[0] -> twin -> prev = hedgehog[n - 1];

		report (hedgehog[n - 1]);
		report (hedgehog[0] -> twin);
		cout << endl;
	}

	set <halfEdge*> went;

	vector <Edge> edgeClosest;
	vector <Pt> tryClosest;
	vector < pair <Pt, Face*> > waitClosest;

	for (halfEdge* a : D.e) {
		if (!a -> incidentFace) {
			Face* f = new Face ();
			a -> incidentFace = f;
			D.f.push_back (f);
			f -> leftmost = a;

			halfEdge* x = a -> next;
			cout << "Face: " << a -> incidentFace << " ";
			report (a);

			while (x != a) {
				cout << "Face: " << a -> incidentFace << " ";
				report (x);
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
				while (x != a) {
					if (x -> p < x -> q) edgeClosest.push_back (Edge (a -> p, a -> q, a));
					x = x -> next;
				}

				tryClosest.push_back (f -> leftmost -> origin -> p);
				waitClosest.push_back ({f -> leftmost -> origin -> p, f});
			}
			else {
				if (a -> p > a -> q) edgeClosest.push_back (Edge (a -> q, a -> p, a));
				x = a -> next;
				while (x != a) {
					if (x -> p > x -> q) edgeClosest.push_back (Edge (a -> q, a -> p, a));
					x = x -> next;
				}
			}
		}
	}

	return D;
}