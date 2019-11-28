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
}
