
#include "lane/polyfit.hpp"
#include <cmath>

namespace lane
{

	double polySingleVal(const vector<double>& p_seq, const double x)
	{
		double y = 0;
		for(size_t i = 0; i < p_seq.size(); i++)
		{
			if(x != 0)
				y += p_seq[i] * GM_Pow(x, i);
		}
		return y;
	}

	vector<double> polyval(const vector<double>& p_seq, const vector<double>& x_seq)
	{
		vector<double> z_seq(x_seq.size());
		for(int n=0; n<x_seq.size(); n++)
		{
			double x = x_seq[n];
			int i;
			double y = 0;
			for(size_t i = 0; i < p_seq.size(); i++)
			{
				if(x != 0)
					y += p_seq[i] * GM_Pow(x, i);
			}
			z_seq[n]=y;
		}
		return z_seq;
	}

	double polyerr(const vector<double> &x_seq, const vector<double> &y_seq, vector<double>& p_seq)	   {
		vector<double> y_val = polyval(p_seq, x_seq);
		double mean_err = 0;
		for (int i = 0; i < x_seq.size(); i++) {
			mean_err += fabs(y_val[i] - y_seq[i]);
		}    
		return mean_err / x_seq.size();
	}

	vector<double> polyfit(const vector<double> &x_seq, const vector<double> &y_seq, int nDegree)
	{
		vector<double> p_seq;
		if(x_seq.size()<2)
		{
			return p_seq;
		}
		if(nDegree > x_seq.size()-1)
		{
			nDegree = x_seq.size()-1;
		}
		if(x_seq.size()!=y_seq.size())
		{
			cerr<<"Error: Dimension of x and y mis match "<<endl;
			return p_seq;
		}
		p_seq.resize(nDegree+1);

		const int nPoints = x_seq.size();
		GMtype_Data DataX, DataY;
		DataX.size = nPoints;
		DataY.size = nPoints;

		for(size_t i=0; i<nPoints; i++)
		{
			DataX.element[i] = x_seq[i];
			DataY.element[i] = y_seq[i];
		}

		GMtype_Polynomial Polynomial;
		Polynomial.degree = nDegree;

		GM_PolyFit(DataX, DataY, &Polynomial);

		for(size_t i=0; i<p_seq.size(); i++)
		{
			p_seq[i]=Polynomial.coef[i];
		}

		return p_seq;
	}

};
