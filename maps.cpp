#include "maps.h"

int main () {
	int n;
	cin >> n;

	vector <Edge> v;

	for (int i = 0; i < n; i++) {
		double x1, y1, x2, y2;
		cin >> x1 >> y1 >> x2 >> y2;
		v.push_back (Edge (x1, y1, x2, y2));
	}

	auto A = construct (v);

	for (Face* f : A.f) {
		for (auto k : f -> inner) {
			cout << f << ' ';
			report (k);
			for (halfEdge* e = k -> next; e != k; e = e -> next) {
				cout << f << ' ';
				report (e);
			}
		}
		cout << f << ' ';
		report (f -> outer);
		for (halfEdge* e = f -> outer -> next; e != f -> outer; e = e -> next) {
			cout << f << ' ';
			report (e);
		}
	}
}
