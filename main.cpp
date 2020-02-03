#include "DCELHandler.h"

int main () {

	DCELHandler make;

	int n, m;
	cin >> n;

	vector <Edge> v;
	vector <Pt> p;

	for (int i = 0; i < n; i++) {
		double x1, y1, x2, y2;
		cin >> x1 >> y1 >> x2 >> y2;
		v.push_back (Edge (x1, y1, x2, y2));
	}

	cin >> m;

	for (int i = 0; i < m; i++) {
		double x, y;
		cin >> x >> y;
		p.push_back (Pt (x, y));
	}

	auto A = make.construct (v);
	make.fill (A, p);

	cin >> n;

	v.clear ();

	for (int i = 0; i < n; i++) {
		double x1, y1, x2, y2;
		cin >> x1 >> y1 >> x2 >> y2;
		v.push_back (Edge (x1, y1, x2, y2));
	}

	cin >> m;
	p.clear ();
	for (int i = 0; i < m; i++) {
		double x, y;
		cin >> x >> y;
		p.push_back (Pt (x, y));
	}

	auto B = make.construct (v);

	make.fill (B, p);

	auto AB = make.merge (A, B);
	cout << "START AB---------------------------------------\n";
	cout << A.outer_bound << ' ' << B.outer_bound << endl;

	for (Face* f : A.f) {
		if (f -> painted) cout << "Painted!\n";
		for (auto k : f -> inner) {
			cout << "Inner:\n";
			cout << f << ' ';
			make.report (k);
			for (halfEdge* e = k -> next; e != k; e = e -> next) {
				cout << f << ' ';
				make.report (e);
			}
		}
		cout << "Outer:\n";
		cout << f << ' ';
		make.report (f -> outer);
		for (halfEdge* e = f -> outer -> next; e != f -> outer; e = e -> next) {
			cout << f << ' ';
			make.report (e);
		}
	}

	cout << "----------------------------\n";

	for (Face* f : B.f) {
		if (f -> painted) cout << "Painted!\n";
		for (auto k : f -> inner) {
			cout << "Inner:\n";
			cout << f << ' ';
			make.report (k);
			for (halfEdge* e = k -> next; e != k; e = e -> next) {
				cout << f << ' ';
				make.report (e);
			}
		}
		cout << "Outer:\n";
		cout << f << ' ';
		make.report (f -> outer);
		for (halfEdge* e = f -> outer -> next; e != f -> outer; e = e -> next) {
			cout << f << ' ';
			make.report (e);
		}
	}

	cout << "----------------------------\n";

	for (Face* f : AB.f) {
		if (f -> painted) cout << "Painted!\n";
		for (auto k : f -> inner) {
			cout << "Inner:\n";
			cout << f << ' ';
			make.report (k);
			for (halfEdge* e = k -> next; e != k; e = e -> next) {
				cout << f << ' ';
				make.report (e);
			}
		}
		cout << "Outer:\n";
		cout << f << ' ';
		make.report (f -> outer);
		for (halfEdge* e = f -> outer -> next; e != f -> outer; e = e -> next) {
			cout << f << ' ';
			make.report (e);
		}
	}
}
