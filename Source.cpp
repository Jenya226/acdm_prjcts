#include <iostream>
#include "kernel.h"
using namespace std;
void menu(int nodes, int kernels, kernel* arr)
{
	double* arr_nodes_X = new double[nodes + 1]{ 0 };
	double* arr_nodes_Y = new double[nodes + 1]{ 0 };
	for (int i = 0; i < nodes; i++)
	{
		cout << " X coordinate of " << i + 1 << " node : ";
		cin >> arr_nodes_X[i];
		cout << " Y coordinate of " << i + 1 << " node : ";
		cin >> arr_nodes_Y[i];
		cout << endl;
	}
	for (int i = 0; i < kernels; i++)
	{
		int temp1, temp2;
		cout << " Enter connections for " << i + 1 << " kernel " << endl;
		cout << " From ";
		cin >> temp1;
		cout << " to ";
		cin >> temp2;
		if (temp1 > temp2)
			swap(temp1, temp2);
		arr[i].link_i = temp1; arr[i].ix = arr_nodes_X[temp1 - 1]; arr[i].iy = arr_nodes_Y[temp1 - 1];
		arr[i].link_j = temp2; arr[i].jx = arr_nodes_X[temp2 - 1]; arr[i].jy = arr_nodes_Y[temp2 - 1];
	}
}
void local_rigidity(double EA, double** arr, kernel* a, int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i].length = sqrt(pow(a[i].jx - a[i].ix, 2) + pow(a[i].jy - a[i].iy, 2));
		a[i].angle = atan((a[i].jy - a[i].iy) / (a[i].jx - a[i].ix));
		double C = EA / a[i].length;
		double local[4][4];
		local[0][0] = local[2][2] = C * cos(a[i].angle) * cos(a[i].angle);
		local[0][2] = local[2][0] = -1 * local[0][0];
		local[0][1] = local[1][0] = local[2][3] = local[3][2] = (C / 2) * sin(2 * a[i].angle);
		local[0][3] = local[3][0] = local[1][2] = local[2][1] = -1 * local[0][1];
		local[1][1] = local[3][3] = C * sin(a[i].angle) * sin(a[i].angle);
		local[3][1] = local[1][3] = -1 * local[1][1];
		int links[4]{ (a[i].link_i * 2 - 1),a[i].link_i * 2,(a[i].link_j * 2 - 1),a[i].link_j * 2 };
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				arr[links[i]][links[j]] = local[i][j] + arr[links[i]][links[j]];


	}
}

double* gauss(double** a, double* y, int n)
{
	double* x, max;
	int k, index;
	const double eps = 0.00001;
	x = new double[n];
	k = 0;
	while (k < n)
	{

		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		if (max < eps)
		{
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue;
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue;
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}
void nulls(double** arr, int a, int b, int n)
{
	arr[0][0] = 1;
	int f;
	int nodes[4]{ a * 2 - 1,a * 2,b * 2 - 1,b * 2 };
	for (int j = 0; j < 4; j++)
	{
		f = nodes[j];
		for (int i = 0; i < n; i++)
		{
			arr[f][i] = arr[i][f] = 0;
			if (i == f)
				arr[i][i] = 1;
		}
	}
}

int main()
{
	double EA = 2 * pow(10, 11) * 0.001;
	int number_of_kernels, number_of_nodes;
	cout << " Quantity of nodes : ";
	cin >> number_of_nodes;
	double** rigidity = new double* [number_of_nodes * 2 + 1];
	for (int i = 0; i < number_of_nodes * 2 + 1; i++)
		rigidity[i] = new double[number_of_nodes * 2 + 1]{ 0 };
	cout << endl; cout << " Quantity of kernels ";
	cin >> number_of_kernels;
	kernel* kernel_array = new kernel[number_of_kernels];
	menu(number_of_nodes, number_of_kernels, kernel_array);
	local_rigidity(EA, rigidity, kernel_array, number_of_kernels);
	for (int i = 0; i < number_of_kernels; i++)
		cout << " Angle :" << kernel_array[i].angle * 57.295779 << " \t length " << kernel_array[i].length << endl;
	int a, b, n; int force_n;
	cout << "Node numbers with anchorage \n First node: ";
	cin >> a;
	cout << "Second node: ";
	cin >> b;
	nulls(rigidity, a, b, number_of_nodes * 2 + 1);
	cout << "Quantity of forces :";
	cin >> force_n;
	double* forces = new double[number_of_nodes * 2 + 1]{ 0 };
	for (int i = 0; i < force_n; i++)
	{
		cout << "Node number : ";
		cin >> n;
		cout << " X: ";
		cin >> forces[n * 2 - 1];
		cout << " Y : ";
		cin >> forces[n * 2];
	}
	forces[a * 2 - 1] = forces[a * 2] = forces[b * 2 - 1] = forces[b * 2 - 1] = 0;
	double* x = gauss(rigidity, forces, number_of_nodes * 2 + 1);
	for (int i = 0; i < number_of_kernels; i++)
	{
		kernel_array[i].deltaL = (x[kernel_array[i].link_j * 2 - 1] - x[kernel_array[i].link_i * 2 - 1]) * cos(kernel_array[i].angle) + (x[kernel_array[i].link_j * 2] - x[kernel_array[i].link_i * 2]) * sin(kernel_array[i].angle);
		cout << "Sigma  " << i + 1 << " = " << (kernel_array[i].deltaL / kernel_array[i].length) * 2 * pow(10, 11) << endl;

	}
	system("pause"); return 0;
}
