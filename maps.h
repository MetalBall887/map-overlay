#include <iostream>
#include "line_segment.cpp"

using namespace std;

DCEL construct (vector <Edge> e) {
	auto res = lineSegInt (e);

	DCEL D;

	map <Pt, vector <halfEdge*> > origins;

	for (auto a : res.second) {
		halfEdge fst (a.p, a.q), snd (a.q, a.p);
		fst.twin = *snd;
		snd.twin = *fst;

		fst.origin = a.q;
		snd.origin = a.p;

		origins[a.q].push_back (fst);
		origins[a.p].push_back (snd);

		D.e.push_back (*fst);
		D.e.push_back (*snd);
	}

	for (auto a : res.first) {
		Vertex vk (a);
		for (auto b : origins[a]) {
			vk.push_back (b);
		}

		auto& hedgehog = vk.incident;
		int n = hedgehog.size ();

		sort (hedgehog.begin (), hedgehog.end (), radial);

		for (int i = 0; i < n - 1; i++) {
			hedgehog[i] -> next = hedgehog[i + 1] -> twin;
			hedgehog[i + 1] -> twin -> prev = hedgehog[i];
		}
		hedgehog[n - 1] -> next = hedgehog[0] -> twin ;
		hedgehog[0] -> twin  -> prev = hedgehog[n - 1];
	}

	set <halfEdge*> went;

	for (halfEdge* a : D.e) {
		if (!went.count (a)) {
			went.insert (a);
			Face* f = new Face ();
			a -> incidentFace = f;
			D.f.push_back (*f);
			f -> leftmost = a -> origin;

			auto x = a -> next;

			while (x != a) {
				x -> incidentFace = a -> incidentFace;
				x = x -> next;
				if (tie (f -> leftmost -> origin -> x, f -> leftmost -> origin -> y) 
					> tie (x -> origin -> x, x -> origin -> y)) {
					f -> leftmost = x;
				}
			}

			if (f -> leftmost -> next == f -> leftmost -> twin || ) {//directed area: bound is outer
				tryClosest.push_back (f -> leftmost -> origin);
				waitClosest[f] = f -> leftmost -> origin;
			} 
		}

		if (a -> p < a -> q) edgeClosest.push_back (Edge (a -> q, a -> q, a -> incidentFace));
	}

	map <Pt, halfEdge>
}