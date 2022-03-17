#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>

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

	//After this while loop, we have all the zero-crossings for negative to positive written in "array_neg"

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
