#include "label_generator.hpp"

void LabelGenerator::outputLabels(ofstream &output_ofs)
{
	int n=0;
	for(n=0; n<labels.size(); n++)
	{
		if(axis_per_line!=labels[n].size())
		{
			std::cerr<<labels[n].size()<<endl;
			std::cerr<<"label generate error"<<endl;
		}
		for(int i=0; i<labels[n].size(); i++)
		{
			output_ofs<<labels[n][i]<<" ";
		}
	}
	// if less than l
	while(n<max_lines)
	{
		for(int i=0; i<axis_per_line; i++)
		{
			output_ofs<<"-1 ";
		}
		n++;
	}
}

void LabelGenerator::outputWeights(ofstream &output_ofs)
{
	int n=0;
	for(n=0; n<weights.size(); n++)
	{
		if(axis_per_line!=weights[n].size())
		{
			std::cerr<<"label generate error"<<endl;
		}
		for(int i=0; i<weights[n].size(); i++)
		{
			output_ofs<<weights[n][i]<<" ";
		}
	}
	// if less than l
	while(n<max_lines)
	{
		for(int i=0; i<axis_per_line; i++)
		{
			output_ofs<<"0 ";
		}
		n++;
	}
}
