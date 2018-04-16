/*************************************************************************
  > File Name: util.cpp
  > Author: Jun Li, Xingang Pan
  > Mail: xingangpan1994@gmail.com 
  > Created Time: 2016年07月20日 星期三 16时20分51秒
 ************************************************************************/

#include "lane/util.hpp"
#include "lane/spline.hpp"

namespace lane
{

	void sortLanes(vector<vector<Point2f> > &lanes)
	{
		for(size_t n=0; n<lanes.size(); n++)
		{
			vector<Point2f> &p_seq = lanes[n];
			sort(p_seq.begin(), p_seq.end(), [](const Point2f& p1, const Point2f &p2){
					return p1.y>p2.y;
					});
			remove_dumplicate_points(p_seq);
		}
		sort(lanes.begin(), lanes.end(), less_by_line_order);
	}

	void nameLineByK(vector<vector<Point2f> > &lanes, vector<int> &lanes_label)
	{	
		if(lanes.empty())
		{
			lanes_label.clear();
			return;
		}
		if(lanes.size()==1)
		{
			if(get_k(lanes[0])<=0)
				lanes_label.push_back(1);
			else
				lanes_label.push_back(2);
		}
		if(lanes.size()==2)
		{
			lanes_label.push_back(1);
			lanes_label.push_back(2);
			return;
		}

		const int max_lines = 4;
		size_t split_idx = lanes.size();
		// lane labels that is less than max_lines are valid
		lanes_label.resize(lanes.size(), max_lines);

		bool find_k_flag = false;
		for(size_t n=1; n<lanes.size(); n++)
		{
			if(get_k(lanes[n-1]) <= 0 && get_k(lanes[n]) > 0)
			{
				split_idx = n;	
				find_k_flag = true;
				break;
			}
		}

		if(!find_k_flag)
		{
			// find from the end_pos
			double center_x = im_width * 0.5;
			for(size_t n=0; n<lanes.size(); n++)
			{
				double interact_x = get_intersection_x(lanes[n]);
				if(interact_x > center_x)
				{
					split_idx = n;	
					break;
				}
			}

		}

		for(int i=0; i<max_lines; i++)
		{
			int idx = i + split_idx - 2;
			if(idx<0 || idx>=lanes.size())
				continue;
			else
				lanes_label[idx] = i;
		}

	}

	void sortLanesCen(vector<vector<Point2f> > &lanes)
	{
		sort(lanes.begin(), lanes.end(), less_by_central_dist);	
		if(lanes.size()<5)
		{	
			sort(lanes.begin(), lanes.end(), less_by_line_order);
			return;
		}
		for(size_t m=lanes.size()-1; m>3; m--)
		{
			lanes.erase(lanes.begin()+m);
		}
		sort(lanes.begin(), lanes.end(), less_by_line_order);
	}
	void remove_dumplicate_points(vector<Point2f> &p_seq)
	{
		int s=0, e=1;
		while(e<p_seq.size())
		{
			if(std::sqrt((p_seq[s].x-p_seq[e].x)*(p_seq[s].x-p_seq[e].x) + (p_seq[s].y-p_seq[e].y)*(p_seq[s].y-p_seq[e].y))<2 || std::fabs(p_seq[s].y-p_seq[e].y) < 1)
			{
				p_seq.erase(p_seq.begin()+e);
			}
			else
			{
				s++;
				e++;
			}
		}
	}

	double get_intersection_x(const vector<Point2f> &line)
	{
		vector<double> x_seq, y_seq;
		for(int i=0; i<line.size(); i++)
		{
			x_seq.push_back(line[i].x);
			y_seq.push_back(line[i].y);
		}
		const int nDegree = 1;
		vector<double> p_seq = polyfit(y_seq, x_seq, nDegree);
		vector<double> y_out(1, im_height);
		vector<double> x_out = polyval(p_seq, y_out);
		return x_out[0];
	}

	double get_k(const vector<Point2f> &line)
	{	
		if(line.size()<2)
		{
			cerr<<"Error: line size must be greater or equal to 2"<<endl;
			return 0;
		}
		return ((line[1].x-line[0].x)/(line[1].y-line[0].y));
	}

	bool less_by_line_order(const vector<Point2f> &line1, const vector<Point2f> &line2)
	{
		double x1 = get_intersection_x(line1);
		double x2 = get_intersection_x(line2);
		return x1<x2;
	}

	bool less_by_central_dist(const vector<Point2f> &line1, const vector<Point2f> &line2)
	{
		double mid = 960;
		double x1 = get_intersection_x(line1);
		double x2 = get_intersection_x(line2);
		return std::abs(x1-mid)<std::abs(x2-mid);
	}


	void changeToFixPoint(vector<Point2f> &curr_lane, size_t fix_pts)
	{
		Spline splineSolver;
		if(curr_lane.size()<2)
			return;
		// change to 30 pts
		vector<Point2f> p_interp = splineSolver.splineInterpStep(curr_lane, 1);		
		if(p_interp.size() < fix_pts)
			return;
		size_t pts_per_line = curr_lane.size();
		double step = (p_interp.size() - 1.0)/(pts_per_line - 1.0);
		curr_lane.clear();
		for(size_t i=0; i<pts_per_line; i++)
		{
			size_t idx = i*step;
			curr_lane.push_back(p_interp[idx]);
		}

	}

	// suppose line has been sorted with y from big to small
	bool extend_line(vector<Point2f> &line)
	{

		double eps = 1;
		if(line.size()<=2)
			return false;
		if(fabs(line[0].x - 0)<eps || fabs(line[0].x - im_width + 1)<eps || fabs(line[0].y - im_height + 1 )<eps)
		{
			return true;
		}
		// roll back first
		int idx = 0;
		while(line[idx].x < 0 || line[idx].x > im_width - 1 || line[idx].y < 0 || line[idx].y > im_height - 1)
		{
			idx++;
			if(idx==line.size())
				break;
		}
		if(idx>line.size()-2)
		{
			cerr<<"error in extending line"<<endl;
			for(auto p:line)
				cerr<<p<<" ";
			cerr<<endl;
			return false;
		}
		Point2f p0 = line[idx];
		Point2f p1 = line[idx+1];
		// check 
		double x0 = p0.x;
		double y0 = p0.y;
		double x1 = p1.x;
		double y1 = p1.y;

		vector<Point2f> new_line;
		if(fabs(x0 - 0)>eps && fabs(x0 - im_width + 1) > eps && fabs(y0 - im_height + 1)>eps)
		{
			double x, y;
			if(fabs(y0 - y1) < 0.001)
			{
				// use x
				if(x0 < x1)
				{
					x = 0;
				}
				else
				{
					x = im_width - 1;
				}
				y = y0 + (x - x0)*(y1-y0)/(x1-x0);
			}
			else
			{
				// first calc y
				y = im_height - 1;
				x = x0 + (y - y0)*(x1 - x0)/(y1-y0); 
				if( x < 0 - eps || x > im_width - 1 + eps)
				{
					// use x
					if(x0 < x1)
					{
						x = 0;
					}
					else
					{
						x = im_width - 1;
					}
					y = y0 + (x - x0)*(y1-y0)/(x1-x0);
				}
			}

			// insert a sequence of point
			Point2f p = Point2f(x,y);
			int insert_num = (std::max(fabs(x - x0), fabs(y - y0)))/eps;
			for(int i=0; i<insert_num; i++)
			{
				double alpha = 1 - i/double(insert_num);
				Point2f tmp_p = alpha*p + (1-alpha)*p0;
				new_line.push_back(tmp_p);
			}
		}

		for(size_t i = idx; i<line.size(); i++)
		{
			new_line.push_back(line[i]);
		}

		std::swap(new_line, line);

		return true;
	}

};
