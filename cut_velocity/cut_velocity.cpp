#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>

// This program writes the averaged cut-velocity signals


//This is the main function for calculating the averaged velocity signals across the zero-crossings
 std::vector<double> cut_velocity(std::vector <double> velocity, std::vector <int> positions, int n)

{
	
	 std::cout << n;  //number of velocity signals each right and left) of the zero crossing

	 std::vector<double> cut_vel;// the array takes in the cut small scale velocity for each zero crossing positions

	int size = positions.size(); // This calculates the size of the positions (of zero-crossing) array, Also needed to calculate the average velocity signal later

	
	for (int i = 0; i < (2*n); i++) //for each velocity signal position
	{

		
		double vel_sum = 0.0; // takes in the sum of the velocity signals for a particular position for all zero crossings

		for (int j = 0; j < size; j++)	//This runs through all the zero-crossings
			
		{
			int temp = positions.at(j); // This gets the positions for each zero crossing window
			int loc = (temp - (n - 1) + i);

			 vel_sum  += velocity.at(loc);
			
		}

		double vel_average = vel_sum /(size * 1.0); // takes the average for each velocity signal position

		cut_vel.push_back(vel_average);

	}


	return cut_vel;

}


int main()
{
	//Read in the positions for positive to negative zero crossings

	std::ifstream infile("pos_neg.out", std::ifstream::binary);

	infile.seekg(0, infile.end);
	int P = infile.tellg();
	infile.seekg(0, infile.beg);

	std::vector<int> array_pos(P / sizeof(int));
	infile.read(reinterpret_cast<char*>(array_pos.data()), static_cast<std::streamsize>(array_pos.size()) * sizeof(double));


	//Read in the positions for negative to positive zero crossings

	std::ifstream infile_a("neg_pos.out", std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int N = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<int> array_neg(N / sizeof(int));
	infile_a.read(reinterpret_cast<char*>(array_neg.data()), static_cast<std::streamsize>(array_neg.size()) * sizeof(int));


	//Read in the the small-scale velocity signals

	std::ifstream infile_b("us_file.out", std::ifstream::binary);

	infile_b.seekg(0, infile_b.end);
	int M = infile_b.tellg();
	infile_b.seekg(0, infile_b.beg);

	std::vector<double> us(M / sizeof(double));
	infile_b.read(reinterpret_cast<char*>(us.data()), static_cast<std::streamsize>(us.size()) * sizeof(double));


	std::vector<double> cut_us_P; // cut small velocity and time for positive to negative
	std::vector<double> cut_us_N; // cut small velocity and time for negative to positive

	int num_of_cuts; // number of velocity signals in the right(or left) of the zero crossing


	// the num_of_cuts variable is read in from a text file
	std::ifstream infile_txt("input.txt");

	infile_txt >> num_of_cuts;

	//call the cut_velocity function to calculate the averaged velocity signals (positive to negative)
	cut_us_P = cut_velocity(us, array_pos, num_of_cuts);

	//call the cut_velocity function to calculate the averaged velocity signals (negative to positive)
	cut_us_N = cut_velocity(us, array_neg, num_of_cuts);

	//write the result to a textfile
	std::ofstream fileOut_P;
	fileOut_P.open("Velocity_pos_neg.txt");

	for (int i = 0; i < cut_us_P.size(); i++)
	{
		fileOut_P << cut_us_P.at(i) << std::endl;
	}
	fileOut_P.close();

	std::ofstream fileOut_N;
	fileOut_N.open("Velocity_neg_pos.txt");

	for (int i = 0; i < cut_us_N.size(); i++)
	{
		fileOut_N << cut_us_N.at(i) << std::endl;
	}
	fileOut_N.close();


	return 0;

}
