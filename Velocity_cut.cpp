

#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>



using Matrix = std::array<std::vector<double>, 100>;
// Allocate enough memory. Will find a way to make this more memory efficient, but works for now

std::vector <Matrix> cut_velocity(std::vector <double> velocity, std::vector <int> positions, int n, std::vector <double> time)

{
	Matrix cut_vel;// the matrix takes in the cut small scale velocity for each zero crossing positions
	Matrix T;// Takes in the time tau for each cut small scale velocity

	int begin = positions.at(0);
	int end = positions.at(positions.size() - 1);

	if (((end + n) > velocity.size()) || ((begin - (n - 1)) < 0)) 
	{
		int n1 = velocity.size() - end;
		int n2 = begin;

		std::cout << "Please choose n that is less than " << (n1) << " and " << (n2) << " " << std::endl;
		exit(-1);
	}

	else
	{
		for (int i = 0; i < positions.size(); i++) // for each of the zero-crossing positions
		{

			int  temp = positions.at(i); // allocates a variable for each zero crossing positions


			for (int j = 0; j < (2 * n); j++) // n is the number of small scale velocity signals and time cut for each zero-crossing position

			{
				int loc = (temp - (n - 1) + j);


				double vel = velocity.at(loc);
				double T_temp = time.at(loc);

				T[i].push_back(T_temp);
				cut_vel[i].push_back(vel);

			}

		}

		return { cut_vel, T };
	}
}


int main()
{

	std::ifstream infile("time_file.out", std::ifstream::binary);

	infile.seekg(0, infile.end);
	int N = infile.tellg();
	infile.seekg(0, infile.beg);

	std::vector<double> t(N / sizeof(double));
	infile.read(reinterpret_cast<char*>(t.data()), t.size() * sizeof(double));


	std::ifstream infile_a("ul_file.out", std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int P = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<double> ul(P / sizeof(double));
	infile_a.read(reinterpret_cast<char*>(ul.data()), ul.size() * sizeof(double));



	std::ifstream infile_b("us_file.out", std::ifstream::binary);

	infile_b.seekg(0, infile_b.end);
	int M = infile_b.tellg();
	infile_b.seekg(0, infile_b.beg);

	std::vector<double> us(M / sizeof(double));
	infile_b.read(reinterpret_cast<char*>(us.data()), us.size() * sizeof(double));

	//std::cout << us.size();



	int i = 1; // to avoid the first "zero-crossing at time = 0"
	std::vector<int> array_pos; //positions for positive to negative

	while (i < ul.size() - 1)
	{
		if (((ul.at(i)) > 1e-12) && ((ul.at(i + 1)) < 1e-12))
		{
			array_pos.push_back(i);

		}

		i++;

	}
	//std::cout << ul.at(45);
	std::cout << " " << std::endl;
	std::cout << " Size of the array for positions with zero-crossing postive to negative is " << array_pos.size() << std::endl;


	int j = 1;
	std::vector<int> array_neg; // to avoid the first "zero-crossing at time = 0"

	while (j < ul.size() - 1)
	{
		if (((ul.at(j)) < 1e-12) && ((ul.at(j + 1)) > 1e-12))
		{
			array_neg.push_back(j);

		}
		j++;
	}

	std::cout << " " << std::endl;
	std::cout << " size of the array for positions with zero-crossing negative to positive is " << array_neg.size() << std::endl;

	//std::cout << array_neg.at(array_neg.size() - 1);
	//std::cout << array_pos.at(array_pos.size() - 1);


	std::vector <Matrix> cut_us_P; // cut small velocity and time for positive to negative
	std::vector <Matrix> cut_us_N; // cut small velocity and time for negative to positive

	int num_of_cuts; // number of velocity signals in the right(or left) of the zero crossing
	std::cout << " " << std::endl;
	std::cout << " How many velocity signals each do you want on the right and left of each zero crossing? ";

		std::cin >> num_of_cuts;


	cut_us_P = cut_velocity(us, array_pos, num_of_cuts, t);

	cut_us_N = cut_velocity(us, array_neg, num_of_cuts, t);

	//std::cout << cut_us_P[1][6][36] << std::endl;

	return 0;



}

