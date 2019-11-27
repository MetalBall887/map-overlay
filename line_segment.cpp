#include "line_segment.h"

int main () {
	int n, m;
	cin >> n >> m;

	vector <Edge> v;
	vector <Pt> p;
	cur = Pt (0, 0);
	cout << (Edge (0, 0, 6, 0) < Edge (3, 0, 3, 5)) << endl;

	for (int i = 0; i < n; i++) {
		double x1, y1;
		cin >> x1 >> y1;
		p.push_back (Pt (x1, y1));
	}

	for (int i = 0; i < m; i++) {
		double x1, y1, x2, y2;

		cin >> x1 >> y1 >> x2 >> y2;

		v.push_back (Edge (x1, y1, x2, y2));
	}

	auto res = lineSegInt (v);

	for (auto a : res.second) {
		cout << a.p.x << ' ' << a.p.y << ' ' << a.q.x << ' ' << a.q.y << endl; 
	}
}
