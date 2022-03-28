#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>
#include <string>
#include <stdio.h>

// This program writes out the positions of zero-crossings in a large-scale velocity signal 

int main()

{

	//This reads in the large scale velocity signals from a binary file and write it into an array
	std::ifstream infile_a("ul_file.out", std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int P = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<double> ul(P / sizeof(double));
	infile_a.read(reinterpret_cast<char*>(ul.data()), static_cast<std::streamsize>(ul.size()) * sizeof(double));


	//This reads in the time signals from a binary file and write it into an array
	std::ifstream infile_t("time_file.out", std::ifstream::binary);

	infile_t.seekg(0, infile_a.end);
	int T = infile_t.tellg();
	infile_t.seekg(0, infile_t.beg);

	std::vector<double> t(T / sizeof(double));
	infile_t.read(reinterpret_cast<char*>(t.data()), static_cast<std::streamsize>(t.size()) * sizeof(double));


	int i = 1; // to avoid the first "zero-crossing at time = 0"
	std::vector<int> array_pos; //positions for positive to negative is written to this array

	while (i < ul.size() - 1)
	{ // 'at' is for assessing an array element 
		if (((ul.at(i)) > 1e-12) && ((ul.at(i + 1)) < 1e-12))
		{
			//for every zero-crossing found, the positions is written to the array "array_pos"
			array_pos.push_back(i);

		}
		i++;

	}

	//After this while loop, we have all the zero-crossings for positive to negative written in "array_pos"

	int j = 1; // to avoid the first "zero-crossing at time = 0"
	std::vector<int> array_neg; //positions for negative to positive is written to this array

	while (j < (ul.size() - 1))
	{
		if (((ul.at(j)) < 1e-12) && ((ul.at(j + 1)) > 1e-12))
		{
			array_neg.push_back(j);

		}
		j++;
	}

	// To define the method for the zero-crossing
	//1 -> Use all zero-crossings
	//2 -> implement a threshold system ( use some zero-crossings)

	int method; // the method you wish to use


	std::ifstream infile_txt("selection.txt"); // read in the method assigned

	infile_txt >> method;

	//If "2" is chosen, the threshold is also read in from a file 
	double threshold;
	std::ifstream infil_txt("criterion.txt");
	infil_txt >> threshold;


	switch (method)
	{

	case 1:

		break;

	case 2:

		std::vector<int> array_neg_;
		std::vector<int> array_pos_;

		array_pos_ = array_pos;
		array_neg_ = array_neg;

		array_pos.clear();
		array_neg.clear();

		for (int i = 1; i < (array_pos_.size()-1); i++)
		{

			double precedent_gap = t[array_pos_.at(i)] - t[array_pos_.at(i - 1)];
			double successive_gap = (t[array_pos_.at(i + 1)] - t[array_pos_.at(i)]);

			if ( ( precedent_gap > (threshold *1.0)) &&  ( successive_gap > (threshold * 1.0)) )
			
			{
				array_pos.push_back(array_pos_.at(i));

			 }


		}

		for (int j = 1; j < (array_neg_.size() -1); j++)
		{
			double precedent_gap = t[array_neg_.at(j)] - t[array_neg_.at(j - 1)];
			double successive_gap = (t[array_neg_.at(j + 1)] - t[array_neg_.at(j)]);

			if ((precedent_gap > (threshold * 1.0)) && (successive_gap > (threshold * 1.0)))
			{
				array_neg.push_back(array_neg_.at(j));

			}

		}
		
	}
	
	for (int i = 0; i < array_pos.size(); i++)
	{

		std::cout << array_pos.at(i) << std::endl;
	}

	for (int i = 0; i < array_neg.size(); i++)
	{

		std::cout << array_neg.at(i) << std::endl;
	}



	// This writes array_pos to a binary file
	std::ofstream out_pos;
	out_pos.open("pos_neg.out");

	if (out_pos)
	{
		out_pos.write(reinterpret_cast<char*>(array_pos.data()), static_cast<std::streamsize>(array_pos.size()) * sizeof(int));
		out_pos.close();

	}

	// This writes array_neg to a binary file
	std::ofstream out_neg;
	out_neg.open("neg_pos.out");

	if (out_neg)
	{
		out_neg.write(reinterpret_cast<char*>(array_neg.data()), static_cast<std::streamsize>(array_neg.size()) * sizeof(int));
		out_neg.close();

	}

	return 0;

}