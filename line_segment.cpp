#include "line_segment.h"

int main () {
	int n, m;
	cin >> m;

	vector <Edge> v;
	vector <Pt> p;
	cur = Pt (0, 0);
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
