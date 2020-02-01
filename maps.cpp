#include "maps.h"

int main () {
	int n, m;
	cin >> n >> m;

	vector <Edge> v;

	for (int i = 0; i < n; i++) {
		double x1, y1, x2, y2;
		cin >> x1 >> y1 >> x2 >> y2;
		v.push_back (Edge (x1, y1, x2, y2));
	}

	v.push_back (Edge (1e5, 1e5, -1e5, 1e5));
	v.push_back (Edge (-1e5, 1e5, -1e5, -1e5));
	v.push_back (Edge (-1e5, -1e5, 1e5, -1e5));
	v.push_back (Edge (1e5, -1e5, 1e5, 1e5));

	auto A = construct (v);
	v.clear ();
	vector <Pt> p;

	for (int i = 0; i < m; i++) {
		double x, y;
		cin >> x >> y;
		p.push_back (Pt (x, y));
	}
	p.clear ();

	fill (A, p);

	cin >> n >> m;

	for (int i = 0; i < n; i++) {
		double x1, y1, x2, y2;
		cin >> x1 >> y1 >> x2 >> y2;
		v.push_back (Edge (x1, y1, x2, y2));
	}

	v.push_back (Edge (1e5, 1e5, -1e5, 1e5));
	v.push_back (Edge (-1e5, 1e5, -1e5, -1e5));
	v.push_back (Edge (-1e5, -1e5, 1e5, -1e5));
	v.push_back (Edge (1e5, -1e5, 1e5, 1e5));

	auto B = construct (v);
	v.clear ();

	for (int i = 0; i < m; i++) {
		double x, y;
		cin >> x >> y;
		p.push_back (Pt (x, y));
	}
	p.clear ();

	fill (B, p);

	cout << "START AB---------------------------------------\n";

	auto AB = merge (A, B);

	for (Face* f : AB.f) {
		if (f -> painted) cout << "Painted!\n";
		for (auto k : f -> inner) {
			cout << "Inner:\n";
			cout << f << ' ';
			report (k);
			for (halfEdge* e = k -> next; e != k; e = e -> next) {
				cout << f << ' ';
				report (e);
			}
		}
		cout << "Outer:\n";
		cout << f << ' ';
		report (f -> outer);
		for (halfEdge* e = f -> outer -> next; e != f -> outer; e = e -> next) {
			cout << f << ' ';
			report (e);
		}
	}
}
