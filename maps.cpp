#include "DCELHandler.h"

using namespace std;

bool DCELHandler::radialComp (halfEdge* a, halfEdge* b) {
	return _radialComp (a -> p - a -> q, b -> p - b -> q);
}

DCEL DCELHandler::construct (std::vector <Edge> e) {

	e.push_back (Edge (1e5, 1e5, -1e5, 1e5));
	e.push_back (Edge (-1e5, 1e5, -1e5, -1e5));
	e.push_back (Edge (-1e5, -1e5, 1e5, -1e5));
	e.push_back (Edge (1e5, -1e5, 1e5, 1e5));
	auto res = lineSegInt (e);

	DCEL D;

	map <Pt, std::vector <halfEdge*> > origins;

	for (auto a : res.second) {
		assert (a.p < a.q);
		halfEdge *fst = new halfEdge (a.p, a.q);
		halfEdge *snd = new halfEdge (a.q, a.p);
		fst -> twin = snd;
		snd -> twin = fst;

		origins[a.q].push_back (fst);
		origins[a.p].push_back (snd);

		D.e.push_back (fst);
		D.e.push_back (snd);
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

		hedgehog[n - 1] -> next = hedgehog[0] -> twin;
		hedgehog[0] -> twin -> prev = hedgehog[n - 1];
	}

	set <halfEdge*> went;

	vector <halfEdge*> edgeClosest;
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
				if (tie (f -> leftmost -> q.x, f -> leftmost -> q.y) > tie (x -> q.x, x -> q.y)) {
					f -> leftmost = x;
				}
				x = x -> next;
			}

			Pt pa = f -> leftmost -> twin -> q;
			Pt pb = f -> leftmost -> q;
			Pt pc = f -> leftmost -> next -> q;


			if ((pa - pb).cross (pc - pb) > 0 || f -> leftmost -> next == f -> leftmost -> twin) {
				if (a -> p > a -> q) edgeClosest.push_back (a);
				x = a -> next;
				while (x != a) {
					if (x -> p > x -> q) edgeClosest.push_back (x);
					x = x -> next;
				}

				tryClosest.push_back (f -> leftmost -> q);
				waitClosest[f -> leftmost -> q] = f;
			}
			else {
				if (a -> p > a -> q) edgeClosest.push_back (a);
				x = a -> next;
				while (x != a) {
					if (x -> p > x -> q) edgeClosest.push_back (x);
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
		} else { D.outer_bound = waitClosest[it.first]; }
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

		if (x != p) delete x;
	}

	D.f.clear ();

	D.f = new_f;



	return D;
}

DCEL::DCEL (const DCEL& other) {
	if (this == &other) return;
	vector <Edge> v;
	for (auto e : other.e) {
		if (e -> p < e -> q)
			v.push_back (Edge (e -> p, e -> q, NULL));
	}

	DCELHandler make;

	auto A = make.construct (v);
	swap (this -> v, A.v);
	swap (this -> e, A.e);
	swap (this -> f, A.f);
	outer_bound = other.outer_bound;
}

void DCELHandler::fill (DCEL& D, vector <Pt> p) {
	vector <halfEdge*> edges;
	for (Face* f : D.f) {
		auto a = f -> outer;

		if (a -> p > a -> q) edges.push_back (a);
		auto x = a -> next;

		while (x != a) {
			if (x -> p > x -> q) edges.push_back (x);
			x = x -> next;
		}

		for (auto a : f -> inner) {
			if (a -> p > a -> q) edges.push_back (a);
			auto x = a -> next;

			while (x != a) {
				if (x -> p > x -> q) edges.push_back (x);
				x = x -> next;
			}
		}
	}

	auto res = findClosest (edges, p);

	for (auto a : res) {
		a.second -> incidentFace -> painted = true;
	}
}

vector <double> DCELHandler::get_line (Pt p, Pt q) {
	double A = p.y - q.y;
	double B = q.x - p.x;
	double C = -p.x * A - p.y * B;

	if (A) B /= A, C /= A, A /= A;
	else C /= B, B /= B;

	return vector <double> ({A, B, C});
}

DCEL DCELHandler::merge (DCEL& A, DCEL& B) {
	vector <Edge> v;

	map <vector <double>, vector <division_data>, comp_triplets> m;
	map <Pt, halfEdge*> from_start;
	vector <Pt> check_a, check_b;
	int cnt = 0;

	halfEdge *is_a = NULL, *is_b = NULL;

	for (auto e : A.e) {
		if (e -> p < e -> q)
			v.push_back (Edge (e -> p, e -> q, NULL));
	}

	for (auto e : B.e) {
		if (e -> p < e -> q)
			v.push_back (Edge (e -> p, e -> q, NULL));
	}

	DCEL AB = construct (v);

	for (auto x : A.e) {
		if (x -> p > x -> q) continue;
		vector <double> v = get_line (x -> p, x -> q);

		m[v].emplace_back (x -> p, 1, x);
		m[v].emplace_back (x -> q, 0, x);
	}

	for (auto x : B.e) {
		if (x -> p > x -> q) continue;
		vector <double> v = get_line (x -> p, x -> q);

		m[v].emplace_back (x -> p, 3, x);
		m[v].emplace_back (x -> q, 2, x);
	}

	for (auto x : AB.e) {
		if (x -> p > x -> q) continue;
		vector <double> v = get_line (x -> p, x -> q);
		m[v].emplace_back ((x -> p + x -> q) / 2.0, 4, x);
	}

	for (auto a : m) {
		auto& v = a.second;

		sort (v.begin(), v.end());

		for (auto x : v) {
			if (x.flag == 0) {
				is_a = NULL;
			}
			if (x.flag == 1){
				is_a = x.ref;
			}
			if (x.flag == 2) {
				is_b = NULL;
			}
			if (x.flag == 3) {
				is_b = x.ref;
			}
			if (x.flag == 4) {
				x.ref -> to_a = (is_a ? is_a -> incidentFace : NULL);
				x.ref -> twin -> to_a = (is_a ? is_a -> twin -> incidentFace : NULL);
				x.ref -> to_b = (is_b ? is_b -> incidentFace : NULL);
				x.ref -> twin -> to_b = (is_b ? is_b -> twin -> incidentFace : NULL);
			}
		}
	}

	map <Pt, Face*> wait_a, wait_b;

	for (auto f : AB.f) {
		halfEdge* l = f -> leftmost;
		f -> to_a = l -> next -> to_a;
		f -> to_b = l -> next -> to_b;
		if (l -> to_a) f -> to_a = l -> to_a;
		if (l -> to_b) f -> to_b = l -> to_b; 

		if (!f -> to_a) {
			wait_a[(l -> p + l -> q) / 2.0] = f;
			check_a.push_back ((l -> p + l -> q) / 2.0);
		}

		if (!f -> to_b) {
			wait_b[(l -> p + l -> q) / 2.0] = f;
			check_b.push_back ((l -> p + l -> q) / 2.0);
		}
	}

	for (int i = 0; i < 2; i++) {
		vector <halfEdge*> edges;
		DCEL* current = (i ? &A : &B);

		for (Face* f : current -> f) {
			auto a = f -> outer;

			if (a -> p > a -> q) edges.push_back (a);
			auto x = a -> next;

			while (x != a) {
				if (x -> p > x -> q) edges.push_back (x);
				x = x -> next;
			}

			for (auto a : f -> inner) {
				if (a -> p > a -> q) edges.push_back (a);
				auto x = a -> next;

				while (x != a) {
					if (x -> p > x -> q) edges.push_back (x);
					x = x -> next;
				}
			}
		}

		auto res = findClosest (edges, (i ? check_a : check_b));

		if (i) {
			for (auto x : wait_a) {
				assert (x.second -> to_a == NULL);
				if (res[x.first]) x.second -> to_a = res[x.first] -> incidentFace;
				else x.second -> to_a = A.outer_bound;
			}
		} else {
			for (auto x : wait_b) {
				assert (x.second -> to_b == NULL);
				if (res[x.first]) x.second -> to_b = res[x.first] -> incidentFace;
				else x.second -> to_b = B.outer_bound;
			}
		}
	}

	for (auto f : AB.f)  {
		f -> painted = (f -> to_a -> painted) && (f -> to_b -> painted);
	}

	return AB;
}