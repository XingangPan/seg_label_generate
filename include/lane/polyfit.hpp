#include <iostream>
#include <vector>

#include "lane/GM_SubLibrary.hpp"
using namespace std;

namespace lane
{

	// take care that we usually set y as x for lanes
	
	double polySingleVal(const vector<double>& p_seq, const double x);

	vector<double> polyval(const vector<double>& p_seq, const vector<double>& x_seq);

	vector<double> polyfit(const vector<double> &x_seq, const vector<double> &y_seq, int nDegree);

	double polyerr(const vector<double> &x_seq, const vector<double> &y_seq, vector<double>& p_seq);

};
