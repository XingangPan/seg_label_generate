#ifndef UTIL_HPP
#define UTIL_HPP

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include "lane/polyfit.hpp"

using namespace std;
using namespace cv;

namespace lane{

const int im_height = 590;
const int im_width = 1640;

void sortLanes(vector<vector<Point2f> > &lanes);
void sortLanesCen(vector<vector<Point2f> > &lanes);
void nameLineByK(vector<vector<Point2f> > &lanes, vector<int> &lanes_label);
void remove_dumplicate_points(vector<Point2f> &p_seq);
// the default value is very important
double get_intersection_x(const vector<Point2f> &line);
bool less_by_line_order(const vector<Point2f> &line1, const vector<Point2f> &line2);
bool less_by_central_dist(const vector<Point2f> &line1, const vector<Point2f> &line2);
double get_k(const vector<Point2f> &line);
void changeToFixPoint(vector<Point2f> &curr_line, size_t fix_pts);
bool extend_line(vector<Point2f> &line);

};

#endif
