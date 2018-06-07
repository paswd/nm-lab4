#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "types.h"

using namespace std;

const string INCORRECT_SELECTION = "Неверная опция";

/*
=====================
FIRST TASK PARAMETERS
=====================
*/

const TNum INTERVAL_BEGIN = 1.;
const TNum INTERVAL_END = 2.;
const TNum INTERVAL_DIFF = INTERVAL_END - INTERVAL_BEGIN;
const TNum y1_first = 1 + exp(.5);
const TNum z1_first = 2 * exp(.5) - 1;
const size_t ORDER = 4;
//const string SPACE = "\t\t";

TNum func_f_first(TNum x, TNum y, TNum z) {
	return ((x*x - 1) * z / x) + ((x*x + 1) * y / (x*x));
}
TNum func_g_first(TNum x, TNum y, TNum z) {
	return 0*x + 0*y + z;
}
TNum func_first_check(TNum x) {
	return (1 + exp(.5*x)) / x;
}

/*
======================
SECOND TASK PARAMETERS
======================
*/

void EulerMethod(void) {
	TNum step;
	cout << "Введите шаг:" << endl;
	cin >> step;
	size_t size = INTERVAL_DIFF / step;
	//cout << fmod(INTERVAL_DIFF, step) << endl;
	cout << "SIZE = " << size << endl;
	vector <TNum> x(size);
	vector <TNum> y(size);
	vector <TNum> z(size);
	vector <TNum> check(size);
	x[0] = INTERVAL_BEGIN;
	y[0] = y1_first;
	z[0] = z1_first;
	check[0] = func_first_check(x[0]);
	size_t runge_size = size / 2;
	vector <TNum> y_runge(runge_size);
	vector <TNum> z_runge(runge_size);
	vector <TNum> runge_res(runge_size);
	y_runge[0] = y[0];
	z_runge[0] = z[0];
	runge_res[0] = 0.;
	cout << "x\t|y\t|z\t|f(xyz)\t|y_r\t|z_r\t|runge\t|y check" << endl; 

	for (size_t i = 1; i < size; i++) {
		x[i] = x[i - 1] + step;
		z[i] = z[i - 1] + step * func_f_first(x[i - 1], y[i - 1], z[i - 1]);
		y[i] = y[i - 1] + step * z[i - 1];
		
		check[i] = func_first_check(x[i]);
	}
	for (size_t i = 1; i < size / 2; i++) {
		y_runge[i] = y_runge[i - 1] + 2 * step * z_runge[i - 1];
		z_runge[i] = z_runge[i - 1] + 2 * step * func_f_first(x[2*i - 2], y_runge[i - 1], z_runge[i - 1]);
		runge_res[i] = (y[2*i] - y_runge[i]) / (pow(2, ORDER) - 1);
	}
	for (size_t i = 0; i < size; i++) {
		cout << x[i] << "\t|";
		cout << y[i] << "\t|";
		cout << z[i] << "\t|";
		cout << func_f_first(x[i], y[i], z[i]) << "\t|";
		if (i % 2 != 0) {
			cout << "---\t|---\t|---\t|";
		} else {
			cout << y_runge[i / 2] << "\t|";
			cout << z_runge[i / 2] << "\t|";
			cout << runge_res[i / 2] << "\t|";
		}
		cout << check[i] << endl;
	}
}

void RungeKuttMethod(void) {
	TNum step;
	cout << "Введите шаг:" << endl;
	cin >> step;
	size_t size = INTERVAL_DIFF / step;
	//cout << fmod(INTERVAL_DIFF, step) << endl;
	cout << "SIZE = " << size << endl;
	vector <TNum> x(size);
	vector <TNum> y(size);
	vector <TNum> z(size);
	vector <TNum> dy(size);
	vector <TNum> dz(size);
	vector <TNum> check(size);
	vector <TNum[4]> k(size);
	vector <TNum[4]> l(size);
	vector <TNum[2]> theta(size);

	x[0] = INTERVAL_BEGIN;
	y[0] = y1_first;
	z[0] = z1_first;

	for (size_t i = 0; i < size - 1; i++) {
		k[i][0] = step * func_f_first(x[i], y[i], z[i]);
		l[i][0] = step * func_g_first(x[i], y[i], z[i]);
		k[i][1] = step * func_f_first(x[i] + .5 * step, y[i] + .5 * k[i][0], z[i] + .5 * l[i][0]);
		l[i][1] = step * func_g_first(x[i] + .5 * step, y[i] + .5 * k[i][0], z[i] + .5 * l[i][0]);
		k[i][2] = step * func_f_first(x[i] + .5 * step, y[i] + .5 * k[i][1], z[i] + .5 * l[i][1]);
		l[i][2] = step * func_g_first(x[i] + .5 * step, y[i] + .5 * k[i][1], z[i] + .5 * l[i][1]);
		k[i][3] = step * func_f_first(x[i] + step, y[i] + k[i][2], z[i] + k[i][2]);
		l[i][3] = step * func_g_first(x[i] + step, y[i] + k[i][2], z[i] + k[i][2]);
		theta[i][0] = abs((k[i][1] - k[i][2]) / (k[i][0] - k[i][1]));
		theta[i][1] = abs((l[i][1] - l[i][2]) / (l[i][0] - l[i][1]));

		dy[i] = (k[i][0] + 2*k[i][1] + 2*k[i][2] + k[i][3]) / 6.;
		dz[i] = (l[i][0] + 2*l[i][1] + 2*l[i][2] + l[i][3]) / 6.;
		x[i + 1] = x[i] + step;
		y[i + 1] = y[i] + dy[i];
		z[i + 1] = z[i] + dz[i];
		check[i] = func_first_check(x[i]);
	}
	cout << "x\t|y\t|z\t|theta\t|check" << endl;
	for (size_t i = 0; i < size - 1; i++) {
		cout << x[i] << "\t|";
		cout << y[i] << "\t|";
		cout << z[i] << "\t|";
		cout << theta[i][0] << "\t|";
		cout << check[i] << endl;
	}

}

void AdamsMethod(void) {
	TNum step;
	cout << "Введите шаг:" << endl;
	cin >> step;
	size_t size = INTERVAL_DIFF / step;
	//cout << fmod(INTERVAL_DIFF, step) << endl;
	cout << "SIZE = " << size << endl;
	vector <TNum> x(size);
	vector <TNum> y(size);
	vector <TNum> z(size);
	vector <TNum> dy(size);
	vector <TNum> dz(size);
	vector <TNum> check(size);
	vector <TNum[4]> k(size);
	vector <TNum[4]> l(size);
	vector <TNum[2]> theta(size);

	x[0] = INTERVAL_BEGIN;
	y[0] = y1_first;
	z[0] = z1_first;

	size_t i;
	for (i = 0; i < min((size_t) 3, size - 1); i++) {
		k[i][0] = step * func_f_first(x[i], y[i], z[i]);
		l[i][0] = step * func_g_first(x[i], y[i], z[i]);
		k[i][1] = step * func_f_first(x[i] + .5 * step, y[i] + .5 * k[i][0], z[i] + .5 * l[i][0]);
		l[i][1] = step * func_g_first(x[i] + .5 * step, y[i] + .5 * k[i][0], z[i] + .5 * l[i][0]);
		k[i][2] = step * func_f_first(x[i] + .5 * step, y[i] + .5 * k[i][1], z[i] + .5 * l[i][1]);
		l[i][2] = step * func_g_first(x[i] + .5 * step, y[i] + .5 * k[i][1], z[i] + .5 * l[i][1]);
		k[i][3] = step * func_f_first(x[i] + step, y[i] + k[i][2], z[i] + k[i][2]);
		l[i][3] = step * func_g_first(x[i] + step, y[i] + k[i][2], z[i] + k[i][2]);
		theta[i][0] = abs((k[i][1] - k[i][2]) / (k[i][0] - k[i][1]));
		theta[i][1] = abs((l[i][1] - l[i][2]) / (l[i][0] - l[i][1]));

		dy[i] = (k[i][0] + 2*k[i][1] + 2*k[i][2] + k[i][3]) / 6.;
		dz[i] = (l[i][0] + 2*l[i][1] + 2*l[i][2] + l[i][3]) / 6.;
		x[i + 1] = x[i] + step;
		y[i + 1] = y[i] + dy[i];
		z[i + 1] = z[i] + dz[i];
		check[i] = func_first_check(x[i]);
	}
	for (; i < size - 1; i++) {
		x[i + 1] = x[i] + step;
		y[i + 1] = y[i] + (step *
			(
				55. * func_g_first(x[i], y[i], z[i])
				- 59. * func_g_first(x[i - 1], y[i - 1], z[i - 1])
				+ 37. * func_g_first(x[i - 2], y[i - 2], z[i - 2])
				- 9. * func_g_first(x[i - 3], y[i - 3], z[i - 3])
			)) / 24.;
		z[i + 1] = z[i] + (step *
			(
				55. * func_f_first(x[i], y[i], z[i])
				- 59. * func_f_first(x[i - 1], y[i - 1], z[i - 1])
				+ 37. * func_f_first(x[i - 2], y[i - 2], z[i - 2])
				- 9. * func_f_first(x[i - 3], y[i - 3], z[i - 3])
			)) / 24.;
	}
	cout << "x\t|y\t|z\t|check" << endl;
	for (i = 0; i < size - 1; i++) {
		check[i] = func_first_check(x[i]);
		cout << x[i] << "\t|";
		cout << y[i] << "\t|";
		cout << z[i] << "\t|";
		cout << check[i] << endl;
	}
}




/*
====
MAIN
====
*/

int main(void) {
	cout.precision(2);
	cout << "=================" << endl;
	cout << "Выберите задание:" << endl;
	cout << "1 - Метод Эйлера" << endl;
	cout << "2 - Метод Рунге-Кутты" << endl;
	cout << "3 - Метод Адамса" << endl;
	cout << "4 - Метод стрельбы" << endl;
	cout << "5 - Конечно-разностный метод" << endl;
	cout << "=================" << endl;

	cout << "Ваш выбор: ";
	size_t selection = 0;
	//size_t in_selection = 0;
	cin >> selection;
	cout << "=================" << endl;

	/*cout << "Имя файла входных данных:" << endl;
	string filename = "";
	cin >> filename;*/

	switch (selection) {
		case 1:
			EulerMethod();
			break;
		case 2:
			RungeKuttMethod();
			break;
		case 3:
			AdamsMethod();
			//RunSolve();
			break;
		case 4:
			//RunSolve();
			break;
		case 5:
			//RunSolve();
			break;
		case '\0':
			break;
		default:
			cout << INCORRECT_SELECTION << endl;
			break;
	}
	return 0;
}
