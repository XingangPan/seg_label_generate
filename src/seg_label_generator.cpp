/*************************************************************************
  > File Name: seg_label_generator.cpp
  > Author: Jun Li, Xingang Pan
  > Mail: xingangpan1994@gmail.com
  > Created Time: 2016年07月27日 星期三 19时51分43秒
 ************************************************************************/


#include <sstream>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include "seg_label_generator.hpp"
#include "lane/polyfit.hpp"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

void SegLabelGenerator::readLabelFile(const string &file, const string &turn_type_file)
{
	ifstream ifs_label(file, ios::in);
	lanes.clear();
	valid_idx_y_start.clear();
	valid_idx_y_end.clear();
	valid_idx_y_start.resize(20, pts_per_line);
	valid_idx_y_end.resize(20, pts_per_line);
	labels.clear();
	weights.clear();
	labels.resize(20, vector<double>(pts_per_line, -1000));
	weights.resize(20, vector<double>(pts_per_line, 0));
	if(ifs_label.fail())
	{
		return;
	}
	string str_lane;
	while(getline(ifs_label, str_lane))
	{
		stringstream ss;
		ss<<str_lane;
		double t_x, t_y;
		vector<Point2f> curr_lane;
		while(ss>>t_x>>t_y)
		{
			curr_lane.push_back(Point2f(t_x, t_y));
		}

		lane::remove_dumplicate_points(curr_lane);
		std::sort(curr_lane.begin(), curr_lane.end(), [](const Point2f& p1, const Point2f& p2)
				{
				return p1.y>p2.y;
				}
			);
		if(curr_lane.size()<2)
			continue;
		// change to 35 pts
		vector<Point2f> p_interp = splineSolver.splineInterpStep(curr_lane, 1);		
		if(p_interp.size() < pts_per_line)
			continue;
		double step = (p_interp.size() - 1.0)/(pts_per_line - 1.0);
		curr_lane.clear();
		for(size_t i=0; i<pts_per_line; i++)
		{
			size_t idx = i*step;
			curr_lane.push_back(p_interp[idx]);
		}
		lanes.push_back(curr_lane);
	}

	sort(lanes.begin(), lanes.end(), lane::less_by_line_order);

	this->updateLabels(lanes);	

	ifs_label.close();
}

void SegLabelGenerator::remove_double_line(vector<vector<Point2f> > &lanes, double dis_thre)
{
	for (size_t n=0; n<lanes.size(); n++)
	{
		vector<Point2f> curr_lane = lanes[n];
		vector<Point2f> p_interp = splineSolver.splineInterpStep(curr_lane, 1);
		vector<double> x_lane, y_lane;
		for(size_t i=0; i<curr_lane.size(); i++)
		{
			x_lane.push_back(curr_lane[i].x);
			y_lane.push_back(curr_lane[i].y);
		}
		vector<double> poly_seq = lane::polyfit(y_lane, x_lane, 3);
		size_t idx_y = 0;
		bool out_flag = false;
		while(idx_y < pts_per_line && y_seq[idx_y] - eps > p_interp[0].y)
		{
			if(idx_y < pts_per_line-1 && p_interp[0].y > y_seq[idx_y+1])
			{
				double x = lane::polySingleVal(poly_seq, y_seq[idx_y]);
				labels[n][idx_y] = x;
				weights[n][idx_y] = 1;
				out_flag = true;
			}
			idx_y++;
		}
		if(out_flag)
		{
			valid_idx_y_start[n] = idx_y - 1;
		}
		else
		{
			valid_idx_y_start[n] = idx_y;
		}
		while(idx_y < pts_per_line)
		{
			double y = y_seq[idx_y];
			if(y - eps < p_interp[p_interp.size() - 1].y)
			{
				break;
			}
			// find nearest
			double x = find_nearest_point(p_interp, y);
			if(x > -500)
			{
				labels[n][idx_y] = x;
				weights[n][idx_y] = 1;
			}
			else
			{
				break;
			}

			idx_y++;
		}
		valid_idx_y_end[n] = idx_y;

		// polyfit and use sigmoid to calc weight
		const double center_y = valid_idx_y_end[n] + 2;
		const double sigma  = 2.2;
		while(idx_y < pts_per_line)
		{
			double x = lane::polySingleVal(poly_seq, y_seq[idx_y]);
			double curr_weight = 1 - 1/(1 + exp(-sigma*(idx_y - center_y)));
			labels[n][idx_y] = x;
			weights[n][idx_y] = curr_weight;
			idx_y++;
		}
	}
	double sum = 0;
	double dis = 0;
	double k = 0;
	vector<Point2f> lane_l, lane_r;
	int valid_id_start, valid_id_end, offset=0, numlanes = lanes.size();
	for (size_t n=1; n<numlanes; n++)
	{
		sum = 0;
		valid_id_start = max(valid_idx_y_start[n-1], valid_idx_y_start[n]);
		valid_id_end = min(valid_idx_y_end[n-1], valid_idx_y_end[n]);
		for (size_t i=valid_id_start; i<valid_id_end; i++)
		{
			sum += labels[n][i] - labels[n-1][i];
		}
		dis = sum/(valid_id_end - valid_id_start);
		if (valid_id_start >= valid_id_end)
		{
			dis = 300;
		}
		lane_l = lanes[n-1-offset];
		lane_r = lanes[n-offset];
		k = (lane_l[0].y - lane_l[lane_l.size()-1].y)/(lane_l[0].x - lane_l[lane_l.size()-1].x);
		//cout << "dis: " << dis << " k: " << k << endl;
		if (dis < dis_thre)
		{
			if (k < 0)
			{
				lanes.erase(lanes.begin()+n-1-offset);
				//cout << "erase: " << n-1-offset << endl;
			}
			else
			{
				lanes.erase(lanes.begin()+n-offset);
				//cout << "erase: " << n-offset << endl;
			}
			offset++;
		}
	}
}

void SegLabelGenerator::updateLabels(vector<vector<Point2f> > &lanes)
{
	vector<int> lanes_label;
	this->remove_double_line(lanes, 120);
	valid_idx_y_start.clear();
	valid_idx_y_end.clear();
	valid_idx_y_start.resize(max_lines, pts_per_line);
	valid_idx_y_end.resize(max_lines, pts_per_line);
	labels.clear();
	weights.clear();
	labels.resize(max_lines, vector<double>(pts_per_line, -1000));
	weights.resize(max_lines, vector<double>(pts_per_line, 0));
	lane::nameLineByK(lanes, lanes_label);
	for(size_t l=0; l<lanes.size(); l++)
	{
		int n = lanes_label[l];
		if(n < max_lines)
		{
			const vector<Point2f> &curr_lane = lanes[l];
			// interp
			vector<Point2f> p_interp = splineSolver.splineInterpStep(curr_lane, 1);
			// polyfit
			vector<double> x_lane, y_lane;
			for(size_t i=0; i<curr_lane.size(); i++)
			{
				x_lane.push_back(curr_lane[i].x);
				y_lane.push_back(curr_lane[i].y);
			}
			vector<double> poly_seq = lane::polyfit(y_lane, x_lane, 3);
			
			// find the valid idx_y
			size_t idx_y = 0;
			bool out_flag = false;
			while(idx_y < pts_per_line && y_seq[idx_y] - eps > p_interp[0].y)
			{
				if(idx_y < pts_per_line-1 && p_interp[0].y > y_seq[idx_y+1])
				{
					double x = lane::polySingleVal(poly_seq, y_seq[idx_y]);
					labels[n][idx_y] = x;
					weights[n][idx_y] = 1;
					out_flag = true;
				}
				idx_y++;
			}
			if(out_flag)
			{
				valid_idx_y_start[n] = idx_y - 1;
			}
			else
			{
				valid_idx_y_start[n] = idx_y;
			}


			// set values
			while(idx_y < pts_per_line)
			{
				double y = y_seq[idx_y];
				if(y - eps < p_interp[p_interp.size() - 1].y)
				{
					break;
				}
				// find nearest
				double x = find_nearest_point(p_interp, y);
				if(x > -500)
				{
					labels[n][idx_y] = x;
					weights[n][idx_y] = 1;
				}
				else
				{
					break;
				}

				idx_y++;
			}
			valid_idx_y_end[n] = idx_y;

			// polyfit and use sigmoid to calc weight
			const double center_y = valid_idx_y_end[n] + 2;
			const double sigma  = 2.2;
			while(idx_y < pts_per_line)
			{
				double x = lane::polySingleVal(poly_seq, y_seq[idx_y]);
				double curr_weight = 1 - 1/(1 + exp(-sigma*(idx_y - center_y)));
				labels[n][idx_y] = x;
				weights[n][idx_y] = curr_weight;
				idx_y++;
			}

		}
	
	}
}

void SegLabelGenerator::showLabels(const string &im_name, int width, int wait_time)
{
	Mat im = imread(im_name, 1);
	Mat im2 = imread(im_name, 1);
	const Scalar color_y_seq = Scalar(70, 120, 60);
	const Scalar color_lines[4] = {Scalar(255, 0, 0), Scalar(0, 255, 0), Scalar(255, 0, 255), Scalar(0, 0, 255)};
	const Scalar color_black = Scalar(0, 0, 0);
	const Scalar color_dark = Scalar(139, 139, 0);
	// draw y_seq
	for(size_t l=0; l<pts_per_line; l++)
	{
		line(im, Point(0, y_seq[l]), Point(im.cols-1, y_seq[l]), color_y_seq, 1);
	}
	
	for(int l=0; l<max_lines; l++)
	{
		vector<Point2f> solid_pts;
		vector<Point2f> dash_pts;
		for(size_t i=0; i<pts_per_line; i++)
		{
			// draw scatter points
			if(weights[l][i] > 0.99)
			{
				circle(im, Point2f(labels[l][i], y_seq[i]), 4, color_black, -1);
			}
			else if(weights[l][i]>0.01)
			{
				circle(im, Point2f(labels[l][i], y_seq[i]), 4, color_dark, -1);
			}
			// add points
			if(weights[l][i]>0.5 && i >= valid_idx_y_start[l] && i<valid_idx_y_end[l])
			{
				solid_pts.push_back(Point2f(labels[l][i], y_seq[i]));
			}
			else if(i >= valid_idx_y_end[l])
			{
				dash_pts.push_back(Point2f(labels[l][i], y_seq[i]));
			}
		}

		// only draw solid_pts
		if(solid_pts.size() >= 2)
		{
			vector<Point2f> p_interp_solid = splineSolver.splineInterpStep(solid_pts, 1);
			for(size_t i=0; i<p_interp_solid.size() - 1; i++)
			{
				line(im, p_interp_solid[i], p_interp_solid[i+1], color_lines[l], width);
			}
		}

	}
	resize(im, im, Size(820,295), 0, 0, INTER_NEAREST);

	imshow(window_name, im);
	//namedWindow("origin", 1);
	//imshow("origin", im2);
	waitKey(wait_time);
}

void SegLabelGenerator::outputimLabels(const string &output_path, const string &sub_im_name, int width, bool Flip)
{
	Mat im(im_height, im_width, CV_8UC1, Scalar(0));
	const Scalar color_y_seq = Scalar(70, 120, 60);
	Scalar color_lines[4] = {Scalar(1), Scalar(2), Scalar(3), Scalar(4)};

	const Scalar color_black = Scalar(0, 0, 0);
	const Scalar color_dark = Scalar(139, 139, 0);
	
	for(int l=0; l<max_lines; l++)
	{
		vector<Point2f> solid_pts;
		vector<Point2f> dash_pts;
		for(size_t i=0; i<pts_per_line; i++)
		{
			// add points
			if(weights[l][i]>0.5 && i >= valid_idx_y_start[l] && i<valid_idx_y_end[l])
			{
				solid_pts.push_back(Point2f(labels[l][i], y_seq[i]));
			}
			else if(i >= valid_idx_y_end[l])
			{
				dash_pts.push_back(Point2f(labels[l][i], y_seq[i]));
			}
		}

		// only draw solid_pts
		if(solid_pts.size() >= 2)
		{
			vector<Point2f> p_interp_solid = splineSolver.splineInterpStep(solid_pts, 1);
			for(size_t i=0; i<p_interp_solid.size() - 1; i++)
			{
				if(Flip)
				{
					line(im, p_interp_solid[i], p_interp_solid[i+1], color_lines[3-l], width);
				}
				else
					line(im, p_interp_solid[i], p_interp_solid[i+1], color_lines[l], width);
			}
		}
	
	}
    
	string out_im_name = sub_im_name.substr(0, sub_im_name.find_last_of(".")) + ".png";
	string impath = output_path+out_im_name;
	int path_len = 0;
	struct stat statbuf, buffer;
	int dir_err = 0;
	string path;
	for (string::iterator it=impath.begin()+1; it!=impath.end(); ++it)
	{
		path_len++;
		if(*it=='/')
		{
			path = impath.substr(0, path_len);
			char *cpath = new char [path.length()+1];
			strcpy(cpath, path.c_str());
			if (stat(cpath, &statbuf) != -1)
			{
				if (!S_ISDIR(statbuf.st_mode))
				{
				}
			}
			else
			{
				dir_err = mkdir(cpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				if (dir_err == -1)
				{
					cout << "Error creating directory!" << endl;
					exit(1);
				}
			}
			delete[] cpath;
		}
	}
	Mat imflip;
	if(Flip)
	{
		flip(im, imflip, 1);
		im = imflip;
	}
	if(stat(impath.c_str(), &buffer) != 0)
	{
		imwrite(impath, im);
	//	cout << impath << endl;
	}
}

double SegLabelGenerator::find_nearest_point(const vector<Point2f> &p_interp, const double y)
{
		double x = -1000;
	double min_dis = 1000;
	for(size_t n=0; n<p_interp.size(); n++)
	{
		const Point2f &p = p_interp[n];
		double dis = fabs(p.y - y);
		if(dis < min_dis)
		{
			min_dis = dis;
			x = p.x;
		}
	}
	if(min_dis > eps)
	{
		cerr<<"min_dis: "<<min_dis<<" lager than eps "<<eps<<"!"<<endl;
	}
	return x;
}

void SegLabelGenerator::outputLabels(ofstream &output_ofs)
{
	int n=0;
	for(n=0; n<max_lines;n++)
	{
		if(axis_per_line!=labels[n].size())
		{
			std::cerr<<labels[n].size()<<endl;
			std::cerr<<"Error: label generate error"<<endl;
		}
		if(valid_idx_y_end[n] - valid_idx_y_start[n] < 2 )
		{
			output_ofs<<' '<<0;
		}
		else
		{
			output_ofs<<' '<<1;
		}
	}
}