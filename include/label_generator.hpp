#ifndef LABEL_GENERATOR_HPP
#define LABEL_GENERATOR_HPP

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace std;
using namespace cv;

class LabelGenerator
{
public:
	LabelGenerator()
	{
	}
	virtual void readLabelFile(const string &file, const string &turn_type_file = "")=0;
	void outputLabels(ofstream &output_ofs);
	void outputWeights(ofstream &output_ofs);
protected:
	vector<vector<double> > labels;
	vector<vector<double> > weights;
	int max_lines;
	int axis_per_line;
};

#endif
