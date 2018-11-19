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
const TNum y1_first = 1. + exp(.5);
//const TNum z1_first = 2. * exp(.5) - 1.;
const TNum z1_first = -1.;

const size_t ORDER = 4;
//const string SPACE = "\t\t";

TNum func_f_first(TNum x, TNum y, TNum z) {
	return ((x*x - 1.) * z / x) + ((x*x + 1.) * y / (x*x));
}
TNum func_g_first(TNum x, TNum y, TNum z) {
	return 0*x + 0*y + z;
}
TNum func_first_check(TNum x) {
	return (1. + exp(.5*x*x)) / x;
}

/*
======================
SECOND TASK PARAMETERS
======================
*/

const TNum y1_second = 1. + 4.*log(2.);
const TNum yn_second = -1. + 3.*log(2.);

const TNum ALPHA[2] = {1, 1};
const TNum BETA[2] = {0, 0};

TNum func_f_second(TNum x, TNum y, TNum z) {
	return (2. * y) / (x*x * (x + 1));
}
TNum func_second_check(TNum x) {
	return -1. + (2./x) + ((2. * (x + 1.))/x) * log(abs(x + 1.));
}

TNum pSecond(TNum x) {
    return 2 / x;
}


TNum qSecond(TNum x) {
    return -1;
}


TNum fSecond(TNum x) {
    return 0;
}

void tmaSecond(TVector &res, TVector &a, TVector &b, TVector &c, TVector &d, TNum shape) {
	TVector p(1, -c[0] / b[0]);
    TVector q(1, d[0] / b[0]);
    TVector x(shape + 1, 0);
    for (size_t i = 1; i <= shape; i++) {
        p.push_back(-c[i] / (b[i] + a[i] * p[i - 1]));
        q.push_back((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]));
    }
    for (size_t ii = shape; ii > 0; ii--) {
    	size_t i = ii - 1;
        x[i] = p[i] * x[i + 1] + q[i];
    }
    res.clear();
    for (size_t i = 0; i < x.size() - 1; i++) {
    	res.push_back(x[i]);
    }
}


/*
=======
METHODS
=======
*/

TNum Round(TNum num, size_t signs) {
	size_t tenSt = 1;
	for (size_t i = 0; i < signs; i++) {
		tenSt *= 10;
	}
	long long tmp = (long long) (num * (TNum) tenSt);
	return (TNum) tmp / (TNum) tenSt;
}

void EulerMethod(void) {
	TNum step;
	cout << "Введите шаг:" << endl;
	cin >> step;
	size_t size = INTERVAL_DIFF / step;
	//cout << fmod(INTERVAL_DIFF, step) << endl;
	cout << "SIZE = " << size << endl;
	TVector x(size);
	TVector y(size);
	TVector z(size);
	TVector check(size);
	/*TNum xCurr = INTERVAL_BEGIN;
	TNum yCurr = y1_first;
	TNum zCurr = z1_first;*/
	x[0] = INTERVAL_BEGIN;
	y[0] = y1_first;
	z[0] = z1_first;
	check[0] = func_first_check(x[0]);
	size_t runge_size = size / 2;
	TVector y_runge(runge_size);
	TVector z_runge(runge_size);
	TVector runge_res(runge_size);
	y_runge[0] = y[0];
	z_runge[0] = z[0];
	runge_res[0] = 0.;
	cout << "x\t|y\t|z\t|y_r\t|z_r\t|runge\t|y check" << endl;

	for (size_t i = 1; i < size; i++) {
		x[i] = x[i - 1] + step;
		z[i] = z[i - 1] + step * func_f_first(x[i - 1], y[i - 1], z[i - 1]);
		y[i] = y[i - 1] + step * z[i];
		
		check[i] = func_first_check(x[i]);
	}
	for (size_t i = 1; i < size / 2; i++) {
		y_runge[i] = y_runge[i - 1] + 2 * step * z_runge[i - 1];
		z_runge[i] = z_runge[i - 1] + 2 * step * func_f_first(x[2*i - 2], y_runge[i - 1], z_runge[i - 1]);
		runge_res[i] = (y[2*i] - y_runge[i]) / (pow(2, ORDER) - 1);
	}
	for (size_t i = 0; i < size; i++) {
		cout << Round(x[i], 2) << "\t|";
		cout << Round(y[i], 2) << "\t|";
		cout << Round(z[i], 2) << "\t|";
		//cout << func_f_first(x[i], y[i], z[i]) << "\t|";
		if (i % 2 != 0) {
			cout << "---\t|---\t|---\t|";
		} else {
			cout << Round(y_runge[i / 2], 2) << "\t|";
			cout << Round(z_runge[i / 2], 2) << "\t|";
			cout << Round(runge_res[i / 2], 2) << "\t|";
		}
		cout << check[i] << endl;
	}
	cout << "Epsilon = " << abs(y[size - 1] - check[size - 1]) << endl;
}

void RungeKuttGetSolution(TNum (*func_f)(TNum, TNum, TNum), TNum (*func_g)(TNum, TNum, TNum),
		TVector &x, TVector &y, TVector &z, vector <TNum[2]> &theta,
		TNum start, TNum step, TNum size, TNum yFirst, TNum zFirst) {
	//TNum intervalDiff = end - start;
	//size_t size = intervalDiff / step;
	//x.resize(size);
	//y.resize(size);
	//z.resize(size);
	//theta.resize(size);

	TVector dy(size);
	TVector dz(size);
	//TVector check(size);
	vector <TNum[4]> k(size);
	vector <TNum[4]> l(size);

	x[0] = start;
	y[0] = yFirst;
	z[0] = zFirst;

	//check[0] = func_first_check(x[0]);

	for (size_t i = 0; i < size - 1; i++) {
		l[i][0] = step * func_f(x[i], y[i], z[i]);
		k[i][0] = step * func_g(x[i], y[i], z[i]);
		l[i][1] = step * func_f(x[i] + .5 * step, y[i] + .5 * k[i][0], z[i] + .5 * l[i][0]);
		k[i][1] = step * func_g(x[i] + .5 * step, y[i] + .5 * k[i][0], z[i] + .5 * l[i][0]);
		l[i][2] = step * func_f(x[i] + .5 * step, y[i] + .5 * k[i][1], z[i] + .5 * l[i][1]);
		k[i][2] = step * func_g(x[i] + .5 * step, y[i] + .5 * k[i][1], z[i] + .5 * l[i][1]);
		l[i][3] = step * func_f(x[i] + step, y[i] + k[i][2], z[i] + k[i][2]);
		k[i][3] = step * func_g(x[i] + step, y[i] + k[i][2], z[i] + k[i][2]);
		theta[i][0] = abs((k[i][1] - k[i][2]) / (k[i][0] - k[i][1]));
		theta[i][1] = abs((l[i][1] - l[i][2]) / (l[i][0] - l[i][1]));

		dy[i] = (k[i][0] + 2*k[i][1] + 2*k[i][2] + k[i][3]) / 6.;
		dz[i] = (l[i][0] + 2*l[i][1] + 2*l[i][2] + l[i][3]) / 6.;
		x[i + 1] = x[i] + step;
		y[i + 1] = y[i] + dy[i];
		z[i + 1] = z[i] + dz[i];
		//check[i + 1] = func_first_check(x[i + 1]);
	}
}

void RungeKuttMethod(void) {
	TNum step;
	cout << "Введите шаг:" << endl;
	cin >> step;
	size_t size = INTERVAL_DIFF / step;
	//cout << fmod(INTERVAL_DIFF, step) << endl;
	cout << "SIZE = " << size << endl;
	TVector x(size);
	TVector y(size);
	TVector z(size);
	vector <TNum[2]> theta(size);
	TVector check(size);

	//
	RungeKuttGetSolution(&func_f_first, &func_g_first, x, y, z, theta, 
			INTERVAL_BEGIN, step, size, y1_first, z1_first);

	for (size_t i = 0; i < size; i++) {
		check[i] = func_first_check(x[i]);
	}

	cout << "x\t|y\t|z\t|theta\t|check" << endl;
	for (size_t i = 0; i < size; i++) {
		cout << Round(x[i], 2) << "\t|";
		cout << Round(y[i], 2) << "\t|";
		cout << Round(z[i], 2) << "\t|";
		cout << Round(theta[i][0], 2) << "\t|";
		cout << check[i] << endl;
	}
	cout << "Epsilon = " << abs(y[size - 1] - check[size - 1]) << endl;

}

void AdamsMethod(void) {
	TNum step;
	cout << "Введите шаг:" << endl;
	cin >> step;
	size_t size = INTERVAL_DIFF / step;
	//cout << fmod(INTERVAL_DIFF, step) << endl;
	cout << "SIZE = " << size << endl;
	TVector x(size);
	TVector y(size);
	TVector z(size);
	TVector check(size);
	vector <TNum[2]> theta(size);

	RungeKuttGetSolution(&func_f_first, &func_g_first, x, y, z, theta,
			INTERVAL_BEGIN, step, size, y1_first, z1_first);

	for (size_t i = 3; i < size - 1; i++) {
		//x[i + 1] = x[i] + step;
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
	for (size_t i = 0; i < size; i++) {
		check[i] = func_first_check(x[i]);
	}
	cout << "x\t|y\t|z\t|check" << endl;
	for (size_t i = 0; i < size; i++) {
		check[i] = func_first_check(x[i]);
		cout << Round(x[i], 2) << "\t|";
		cout << Round(y[i], 2) << "\t|";
		cout << Round(z[i], 2) << "\t|";
		cout << check[i] << endl;
	}
	cout << "Epsilon = " << abs(y[size - 1] - check[size - 1]) << endl;
}

TNum derOne(TVector &x, TVector &y, TNum val) {
    size_t i = 0;
    while (x[i + 1] < val - 1e-7 && i < x.size() - 1) {
        i++;
    }
    return (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

TNum next(TVector &x, TVector &yCurr, TVector &yPrev, TNum nCurr, TNum nPrev, TNum alpha, TNum beta) {
	TNum n1 = beta * derOne(x, yCurr, INTERVAL_END);
	TNum n2 = beta * derOne(x, yPrev, INTERVAL_END);
	TNum n3 = alpha * yPrev[x.size() - 1];
	TNum n4 = alpha * yCurr[x.size() - 1] + 1 - yn_second;
	TNum n5 = alpha * yCurr[x.size() - 1] + n1 - n3 - n2;

	return nCurr - n4 * (nCurr - nPrev) / n5;
}

void ShootingMethod() {
	TNum step;
	cout << "Шаг" << endl;
	cin >> step;
	TNum eps;
	cout << "Введите точность" << endl;
	cin >> eps;
	
	//A = y1_second;
	//B = yn_second;

	size_t size = INTERVAL_DIFF / step;

	TNum nPrev = 1.;
	TNum nCurr = .8;

	TVector x(size);
	TVector yPrev(size);
	TVector yCurr(size);
	TVector z(size);
	vector <TNum[2]> theta(size);
	

	RungeKuttGetSolution(&func_f_second, &func_g_first, x, yPrev, z, theta,
			INTERVAL_BEGIN, step, size, nPrev, (y1_second - ALPHA[0] * nPrev) / BETA[0]);
	RungeKuttGetSolution(&func_f_second, &func_g_first, x, yCurr, z, theta,
			INTERVAL_BEGIN, step, size, nCurr, (y1_second - ALPHA[0] * nCurr) / BETA[0]);

	size_t cnt = 0;
	//cout << abs(ALPHA[1] * yCurr[size - 1] +
          //BETA[1] * derOne(x, yCurr, INTERVAL_END) - yn_second) << endl;
	//cout << yCurr[size - 1] << endl
	//cout << ;
	while (abs(ALPHA[1] * yCurr[size - 1] +
          BETA[1] * derOne(x, yCurr, INTERVAL_END) - yn_second) > eps) {
		TNum n = next(x, yCurr, yPrev, nCurr, nPrev, ALPHA[1], BETA[1]);
		nPrev = nCurr;
		nCurr = n;
		yPrev = yCurr;
		RungeKuttGetSolution(&func_f_second, &func_g_first, x, yCurr, z, theta,
			INTERVAL_BEGIN, step, size, nCurr, (y1_second - ALPHA[0] * nCurr) / BETA[0]);
		cnt++;
	}
	cout << cnt << endl;

	cout << "x\t|y\t|z\t|theta\t|check" << endl;
	for (size_t i = 0; i < size; i++) {
		cout << Round(x[i], 2) << "\t|";
		cout << Round(yCurr[i], 2) << "\t|";
		cout << Round(z[i], 2) << "\t|";
		cout << Round(theta[i][0], 2) << "\t|";
		cout << func_second_check(x[i]) << endl;
		//cout << 0 << endl;
	}
	cout << "Epsilon = " << abs(yCurr[size - 1] - func_second_check(x[size - 1])) << endl;

}

void FiniteDifferenceMethod() {
	TVector x(1, INTERVAL_BEGIN);
	TVector a(0);
	TVector b(0);
	TVector c(0);
	TVector d(0);
	TVector y(0);

	TNum step;
	cout << "Шаг" << endl;
	cin >> step;

	size_t size = INTERVAL_DIFF / step;

	a.push_back(0.);
    b.push_back(-2. / (step * (2. - pSecond(INTERVAL_BEGIN) * step)) + qSecond(INTERVAL_BEGIN) * step /
             (2. - pSecond(INTERVAL_BEGIN) * step) + ALPHA[0] / BETA[0]);
    c.push_back(2. / (step * (2. - pSecond(INTERVAL_BEGIN) * step)));
    d.push_back(y1_second / BETA[0] + step * fSecond(INTERVAL_BEGIN) / (2. - pSecond(INTERVAL_BEGIN) * step));
    x.push_back(x[0] + step);

    for (size_t i = 1; i < size; i++) {
    	a.push_back(1. / step*step - pSecond(x[i]) / (2. * step));
        b.push_back(-2. / step*step + qSecond(x[i]));
        c.push_back(1. / step*step + pSecond(x[i]) / (2. * step));
        d.push_back(fSecond(x[i]));
        x.push_back(x[i] + step);
    }

    a.push_back(-2. / (step * (2. + pSecond(x[size - 1]) *  step)));
    b.push_back(2. / (step * (2. + pSecond(x[size - 1]) * step)) - qSecond(x[size - 1]) * step /
             (2. + pSecond(x[size - 1]) * step) + ALPHA[1] / BETA[1]);
    c.push_back(0);
    d.push_back(yn_second / BETA[1] - step * fSecond(x[size - 1]) / (2. + pSecond(x[size - 1]) * step));
    tmaSecond(y, a, b, c, d, a.size());

    cout << "x\t|y\t|check" << endl;
	for (size_t i = 0; i < x.size(); i++) {
		cout << Round(x[i], 2) << "\t|";
		cout << Round(y[i], 2) << "\t|";
		//cout << func_second_check(x[i]) << endl;
		cout << 0 << endl;
	}
	//cout << "Epsilon = " << abs(yCurr[size - 1] - func_second_check(x[size - 1])) << endl;
}




/*
====
MAIN
====
*/

int main(void) {
	//cout.precision(2);
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
			break;
		case 4:
			ShootingMethod();
			break;
		case 5:
			FiniteDifferenceMethod();
			break;
		case '\0':
			break;
		default:
			cout << INCORRECT_SELECTION << endl;
			break;
	}
	return 0;
}
