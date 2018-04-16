/*************************************************************************
	> File Name: main.cpp
	> Author: Jun Li, Xingang Pan
	> Mail: xingangpan1994@gmail.com
	> Created Time: 2016年07月28日 星期四 13时22分56秒
 ************************************************************************/

#include "unistd.h"
#include "seg_label_generator.hpp"
#include <iostream>
using namespace std;

void help()
{
	cout<<"generate seg_label"<<endl;
	cout<<"Usage: "<<endl;
	cout<<"./seg_label_generate [OPTIONS]"<<endl;
	cout<<"-h                 : help"<<endl;
	cout<<"-l                 : list_file"<<endl;
	cout<<"-d                 : dataset path"<<endl;
	cout<<"-s                 : show result"<<endl;
	cout<<"-o                 : output path"<<endl;
	cout<<"-w                 : label width"<<endl;
	cout<<"-i                 : interval"<<endl;
	cout<<"-f                 : start_frame"<<endl;
	cout<<"-p                 : flip"<<endl;
}

int main(int argc, char **argv)
{
	string list_file = "~/works/SCNN/data/CULane/list/train.txt";
	string dir_im = "~/works/SCNN/data/CULane";
	string output_path = "~/works/SCNN/data/CULane/laneseg_label";
	string output_file = "";
	string mode = "imgLabel";   // set mode to "imgLabel" or "trainList"
	bool is_show = false;
	int width = 16;
	int interval = 1;
	int start_frame = 1;
	int set_id = 1;
	// read list
	int ch;
	bool flip = false;
	while((ch = getopt(argc, argv, "hsl:m:d:o:w:i:f:p")) != -1 )
	{
		switch(ch)
		{
			case 'h':
				help();
				return 0;
			case 'l':
				list_file = optarg;
				break;
			case 'm':
				mode = optarg;
				break;
			case 'd':
				dir_im = optarg;
				break;
			case 's':
				is_show = true;
				break;
			case 'o':
				output_path = optarg;
				break;
			case 'w':
				width = atoi(optarg);
				break;
			case 'i':
				interval = atoi(optarg);
				break;
			case 'f':
				start_frame = atoi(optarg);
				break;
			case 'p':
				flip = true;
				break;
		}
	}

	ofstream ofs_out_file;
	if(mode == "trainList")
	{
		output_file = list_file.substr(0, list_file.find_last_of(".")) + "_gt.txt";
		ofs_out_file.open(output_file, ios::out);
	}
	
	SegLabelGenerator seg_label_generator;	

	ifstream ifs_list(list_file, ios::in);
	if(ifs_list.fail())
	{
		cerr<<"file "<<list_file<<" not exist!"<<endl;
		return 1;
	}
	string sub_im_name;
	int count = 0;
	while(getline(ifs_list, sub_im_name))
	{
		count++;
		if(count < start_frame || count%interval != 0)
			continue;
		string im_name = dir_im + sub_im_name;
		string line_label_file = im_name.substr(0, im_name.find_last_of(".")) + ".lines.txt";
		seg_label_generator.readLabelFile(line_label_file);
		if(count%100==0)
			cout << count << ": " << im_name << endl;
		if(is_show)
			seg_label_generator.showLabels(im_name, width);
		else
		{
			// output result
			if(mode == "imgLabel")
			{
				seg_label_generator.outputimLabels(output_path, sub_im_name, width, flip);
			}
			else if(mode == "trainList")
			{
				string segLabel_name = "/laneseg_label_w16/" + sub_im_name.substr(0, sub_im_name.find_last_of(".")) + ".png";
				ofs_out_file<<sub_im_name<<' '<<segLabel_name;
				seg_label_generator.outputLabels(ofs_out_file);
				ofs_out_file<<endl;
			}
			else
			{
				cerr << "illegal mode:" << mode << '!' << endl;
				return 1;
			}
		}
	}

	if(mode == "trainList")
	{
		ofs_out_file.close();
	}
	ifs_list.close();

	return 0;
}
