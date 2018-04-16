#ifndef Y_LABEL_GENERATOR_HPP
#define Y_LABEL_GENERATOR_HPP

#include "label_generator.hpp"
#include "lane/spline.hpp"
#include "lane/util.hpp"

class SegLabelGenerator : public LabelGenerator{
	public:
		SegLabelGenerator()
		{
			im_width = lane::im_width;
			im_height = lane::im_height;
			max_lines = 4;
			// here axis_per_line = 30 because we only output x axis
			axis_per_line = 35;
			pts_per_line = 35;
			y_step = 10;
			eps = 1;
			for(size_t i=0; i<pts_per_line; i++)
			{
				y_seq.push_back(im_height - i*y_step);
			}
			window_name = "im";
			namedWindow(window_name, 1);
		}
		virtual void readLabelFile(const string &file, const string &turn_type_file = "");
		void showLabels(const string &im_name, int width, int wait_time = 0);
		void updateLabels(vector<vector<Point2f> > &lanes);
		void outputLabels(ofstream &output_ofs);
		void outputimLabels(const string &output_path, const string &im_name, int width, bool flip);
		void remove_double_line(vector<vector<Point2f> > &lanes, double dis_thre);

	protected:
		double find_nearest_point(const vector<Point2f> &p_interp, const double y);

	private:
		vector<vector<Point2f> > lanes;
		lane::Spline splineSolver;
		size_t pts_per_line;
		vector<double> y_seq;
		double y_step;
		vector<size_t> valid_idx_y_start;
		vector<size_t> valid_idx_y_end;
		double eps;
		string window_name;
		int im_width;
		int im_height;
};

#endif
