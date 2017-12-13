/*
	File : autotest.cpp
	Aim  : to generate user-specified test suite for programs in Siemens Suite
*/

#include <iostream>
#include <fstream>

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
void random_list(int a[], int n) {
	for (int i = 0; i < n; i++) {
		a[i] = std::rand();
		if (std::rand() % 2)
			a[i] = -a[i];
	}
}
void limits_list(int a[], int n, int min, int max) {
	for (int i = 0; i < n; i++) {
		int val = std::rand();
		val = val % (max - min);
		a[i] = val + min;
	}
}
void output_array(int a[], int n, std::ostream & out) {
	for (int i = 0; i < n; i++) {
		out << a[i] << " ";
	}
	out << std::endl;
}

/* for bubble.c */
/*
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

/* for Day.c */
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

/* for mid.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "mid";
	std::ofstream out(prefix + prname + "/tsuite");

	int a1[3] = { -4, -2, -1 };
	int a2[3] = { 5, 8, 4 };
	int a3[3] = { -5, 0, - 2 };
	int a4[3] = { 1, 2, 0 };
	int a5[3] = { 3, 2, 1 };
	int a6[3] = { 10, 5, 1 };
	int a7[3] = { -5, -1, -3 };
	int a8[3] = {-1, 0, 1};
	int a9[3] = { -1, -5, -3 };
	int a10[3] = { 0, -5, -3 };
	int a11[3] = { -3, -5, 0 };
	int a12[3] = { 8, 9, 3 };
	int a13[3] = { 4, 5, 3 };
	int a14[3] = { -2, 0, 1 };
	int a15[3] = { 5, 7, 0 };
	int a16[3] = { 4, 5, 3 };

	output_permutation(a1, 3, out);
	output_permutation(a2, 3, out);
	output_permutation(a3, 3, out);
	output_permutation(a4, 3, out);
	output_permutation(a5, 3, out);
	output_permutation(a6, 3, out);
	output_permutation(a7, 3, out);
	output_permutation(a8, 3, out);
	output_permutation(a9, 3, out);
	output_permutation(a10, 3, out);
	output_permutation(a11, 3, out);
	output_permutation(a12, 3, out);
	output_permutation(a13, 3, out);
	output_permutation(a14, 3, out);
	output_permutation(a15, 3, out);
	output_permutation(a16, 3, out);

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/* minmax.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "minmax";
	std::ofstream out(prefix + prname + "/tsuite");

	int a1[5] = { 0, 1, 2, 3, 4 };
	int a2[5] = { -5, 7, 19, 80, 888 };
	int a3[5] = { -1 - 1, -1, -2, -3 };
	int a4[5] = { 15, 9, -5, -7, -12 };
	int a5[5] = {0, -1, 1, -1, 0};

	output_permutation(a1, 5, out);
	output_permutation(a2, 5, out);
	output_permutation(a3, 5, out);
	output_permutation(a4, 5, out);
	output_permutation(a5, 5, out);

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/* prime.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "prime";
	std::ofstream out(prefix + prname + "/tsuite");

	int min = -9, max = 727;
	for (int k = 0; k <= max; k++)
		out << k << "\n";

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/* profit.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "profit";
	std::ofstream out(prefix + prname + "/tsuite");

	for (int k = 0; k < 527; k++)
		out << 2000 * k << "\n";

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/* triangle.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "triangle";
	std::ofstream out(prefix + prname + "/tsuite");

	int a1[3] = { 3, 2, 5 };
	int a2[3] = { 5, -1, 3 };
	int a3[3] = { 5, 4, 0 };
	int a4[3] = { 6, 6, 2 };
	int a5[3] = { -9, -10, -7 };
	int a6[3] = { 10, 5, 1 };
	int a7[3] = { -5, -1, -3 };

	output_permutation(a1, 3, out);
	output_permutation(a2, 3, out);
	output_permutation(a3, 3, out);
	output_permutation(a4, 3, out);
	output_permutation(a5, 3, out);
	output_permutation(a6, 3, out);
	output_permutation(a7, 3, out);

	out << std::endl; out.close();
	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/* Calendar.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "Calendar";
	std::ofstream out(prefix + prname + "/tsuite");

	int lower = -405, upper = 871;
	for (int i = lower; i <= upper; i++)
		out << i << "\n";
	out.close();

	return 0;
}
*/

/* new-bubble */
/*
int main() {
	std::string prefix = "../../../MyData/NewBanchmark/";
	//std::string prname = "bubble";
	std::string prname = "insert";
	//std::string prname = "qsort";
	//std::string prname = "minmax";
	std::ofstream out(prefix + prname + "/tsuite");

	const int sizes = 97;
	const int tries = 16;
	const int lists = 32;

	int *list, n, val; 
	for (int i = 0; i < tries; i++) {
		n = std::rand();
		n = std::abs(n);
		n = n % sizes + 1;

		list = new int[n];
		for (int j = 0; j < lists; j++) {
			random_list(list, n);
			output_array(list, n, out);
		}
		
		val = std::rand();
		limits_list(list, n, val, val + 1);
		output_array(list, n, out);

		delete list;
	}

	std::cout << "Press any key to exit...";
	getchar(); exit(0);
}
*/

/* gzip.c */
/*
int main() {
	std::string prefix = "../../../MyData/SiemensSuite/";
	std::string prname = "gzip";
	std::ofstream out(prefix + prname + "/suites/test");

	for (int i = 1; i <= 25; i++) {
		out << "-P[-t]  -I[mydir/file" << i << ".gz]\n";
		out << "-P[-t]  -I[mydir/file" << i << ".z]\n";
		out << "-P[-d]  -I[mydir/file" << i << ".z]\n";
		out << "-P[--decompress]  -I[mydir/file" << i << ".gz]\n";
	}
	
	for (int i = 1; i <= 25; i++) {
		out << "-P[-fqr1]  -I[mydir/file" << i << "]\n";
		out << "-P[-fqr2]  -I[mydir/file" << i << "]\n";
		out << "-P[-fqv3]  -I[mydir/file" << i << "]\n";
		out << "-P[-fqv4]  -I[mydir/file" << i << "]\n";
		out << "-P[-fqrv5]  -I[mydir/file" << i << "]\n";
		out << "-P[-fqrv6]  -I[mydir/file" << i << "]\n";
		out << "-P[-frv7]  -I[mydir/file" << i << "]\n";
		out << "-P[-frv8]  -I[mydir/file" << i << "]\n";
		out << "-P[-fq9]  -I[mydir/file" << i << "]\n";
		out << "-P[-fq1]  -I[mydir/file" << i << "]\n";
		out << "-P[-fr2]  -I[mydir/file" << i << "]\n";
		out << "-P[-fr3]  -I[mydir/file" << i << "]\n";
	}
	out << std::endl; out.close();

	return 0;
}*/