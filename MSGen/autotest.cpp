/*
	File : autotest.cpp
	Aim  : to generate user-specified test suite for programs in Siemens Suite
*/

#include <iostream>
#include <fstream>

/*
void swap_array(int *a, int i, int j) {
	if (i == j) 
		return;
	else {
		int temp;
		temp = a[i];
		a[i] = a[j];
		a[j] = temp;
	}
}
void output_permutation(int a[], int k, 
	int n, std::ostream & out) {
	if (n <= 0) return;

	if (k == 1) {
		for (int i = 0; i < n; i++)
			out << a[i] << " ";
		out << "\n";
	}
	else {
		for (int i = 0; i < k; i++) {
			swap_array(a, n - k, n - k + i);
			output_permutation(a, k - 1, n, out);
			swap_array(a, n - k, n - k + i);
		}
	}
}
void output_permutation(int a[], int n, std::ostream & out) {
	output_permutation(a, n, n, out);
}

int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "bubble";
	std::ofstream out(prefix + prname + "/tsuite");

	int a[5] = { 1, 5, -1, 0, -7 };
	output_permutation(a, 5, out);

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "Day";
	std::ofstream out(prefix + prname + "/tsuite");

	int years[10] = {-1, -5, 0, 400, 100, 700, 32, 128, 79, 21};
	int months[15] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,-6};
	int days[5] = { 1, 0, -1, 27, 16 };

	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 15; j++) {
			for (int k = 0; k < 5; k++) {
				int year = years[i];
				int month = months[j];
				int day = days[k];

				out << year << "\t" << month << "\t" << day << "\n";
			}
		}
	}

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/