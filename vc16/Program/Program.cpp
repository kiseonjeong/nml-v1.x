// Program.cpp : 콘솔 응용 프로그램에 대한 진입점을 정의합니다.
//

#include "stdafx.h"

void test_broadcast();
void test_decomposition();
void test_interpolation();
void test_regression();
void test_spline();

int _tmain(int argc, _TCHAR* argv[])
{
	clock_t t;
	numat m0(msize(4, 4));
	numat m1(msize(4, 4));
	t = clock();
	for (int i = 0; i < 100000; i++)
	{
// 		for (int j = 0; j < m0.rows; j++)
// 		{
// 			for (int k = 0; k < m0.cols; k++)
// 			{
// 				m0(0) = 3.14;
// 			}
// 		}
		numat m2(m0);
	}
	t = clock() - t;
	cout << (float)t << "ms" << endl;
	test_broadcast();
	test_decomposition();
	test_interpolation();
	test_regression();
	test_spline();

	return 0;
}

void test_broadcast()
{
	numat m0(msize(4, 4));
	m0(0, 0) = 100; m0(0, 1) = 3; m0(0, 2) = 45; m0(0, 3) = 33;
	m0(1, 0) = 5; m0(1, 1) = 3; m0(1, 2) = 5; m0(1, 3) = 65;
	m0(2, 0) = 100; m0(2, 1) = 23; m0(2, 2) = 45; m0(2, 3) = 56;
	m0(3, 0) = 5; m0(3, 1) = 3; m0(3, 2) = 9; m0(3, 3) = 8;
	m0.cout();
	m0.row(0, 3, 2).col(0, 3, 2).cout();
	numat m1(msize(4, 1));
	m1(0) = 10;
	m1(1) = 20;
	m1(2) = 30;
	m1(3) = 40;
	m1.cout();
	numat m2 = m0.mul(m1);
	m2.cout();
	numat m3 = m1.mul(m0);
	m3.cout();
	numat m4(msize(1, 4));
	m4(0) = 10;
	m4(1) = 20;
	m4(2) = 30;
	m4(3) = 40;
	m4.cout();
	numat m5 = m0.mul(m4);
	m5.cout();
	numat m6 = m4.mul(m0);
	m6.cout();
	numat m7;
	m7.cout();
}

void test_decomposition()
{
	clock_t t;
	t = clock();

	numat m0(msize(4, 4));
	m0(0, 0) = 100; m0(0, 1) = 3; m0(0, 2) = 45; m0(0, 3) = 33;
	m0(1, 0) = 5; m0(1, 1) = 3; m0(1, 2) = 5; m0(1, 3) = 65;
	m0(2, 0) = 100; m0(2, 1) = 23; m0(2, 2) = 45; m0(2, 3) = 56;
	m0(3, 0) = 5; m0(3, 1) = 3; m0(3, 2) = 9; m0(3, 3) = 8;
	m0.cout();
	LUD lu(m0, LU_PARTIAL_PIVOT);
	lu[0].cout();
	lu[1].cout();
	(lu[0] * lu[1]).cout();
	QRD qr(m0);
	qr[0].cout();
	qr[1].cout();
	(qr[0] * qr[1]).cout();
	(qr[0] * qr[0].T()).cout();
	EVD ev(m0);
	ev[0].cout();
	ev[1].cout();
	(ev[0] * ev[1] * ev[0].I()).cout();
	numat m1(msize(1, 4));
	m1(0) = 5; m1(1) = 3; m1(2) = 9; m1(3) = 8;
	m1.cout();
	SVD sv(m1);
	sv[0].cout();
	sv[1].cout();
	sv[2].cout();
	(sv[0] * sv[1] * sv[2].T()).cout();
	numat m2(msize(3, 3));
	m2(0, 0) = 1; m2(0, 1) = 2; m2(0, 2) = 3;
	m2(1, 0) = 4; m2(1, 1) = 5; m2(1, 2) = 6;
	m2(2, 0) = 7; m2(2, 1) = 8; m2(2, 2) = 9;
	(m2 * SVD::pinv(m2)).cout();
	SVD::pinv(m2).cout();

	t = clock() - t;
	cout << (float)t << "ms" << endl;
}

void test_interpolation()
{
	numat X(msize(6, 1));
	numat Y(msize(6, 1));
	X(0) = 0.0; X(1) = 1.0; X(2) = 2.0; X(3) = 3.0; X(4) = 4.0; X(5) = 5.0;
	Y(0) = 2.1; Y(1) = 7.7; Y(2) = 13.6; Y(3) = 27.2; Y(4) = 40.9; Y(5) = 61.1;
	vandermonde vd(X, Y);
	double result = vd.interpolate(4.5);
	numat C = vd.coeff();
	C.cout();
	newton nt(X, Y);
	result = nt.interpolate(4.5);
	C = nt.coeff();
	C.cout();
	lagrange lg(X, Y);
	result = lg.interpolate(4.5);
}

void test_regression()
{
	clock_t t;
	t = clock();

	numat X(msize(6, 1));
	numat Y(msize(6, 1));
	X(0) = 0.0; X(1) = 1.0; X(2) = 2.0; X(3) = 3.0; X(4) = 4.0; X(5) = 5.0;
	Y(0) = 2.1; Y(1) = 7.7; Y(2) = 13.6; Y(3) = 27.2; Y(4) = 40.9; Y(5) = 61.1;
	curvefit cf(X, Y, 3);
	double result = cf.fit(4.5);
	numat C = cf.coeff();
	C.cout();
	numat Xp = numat::pow(X, 3);
	for (int i = 2; i >= 0; i--)
	{
		Xp = Xp.happend(numat::pow(X, i));
	}
	Xp.cout();
	C = (Xp.T() * Xp).I() * Xp.T() * Y;
	C.cout();

	t = clock() - t;
	cout << (float)t << "ms" << endl;
}

void test_spline()
{
	numat X(msize(4, 1));
	numat Y(msize(4, 1));
	X(0) = 3.0; X(1) = 4.5; X(2) = 7.0; X(3) = 9.0;
	Y(0) = 2.5; Y(1) = 1.0; Y(2) = 2.5; Y(3) = 0.5;
	linear_spline sp1(X, Y);
	double result = sp1.interpolate(5.0);
	numat C = sp1.coeff();
	C.cout();
// 	FILE* writer;
// 	fopen_s(&writer, "D:\\linear_spline.csv", "w");
// 	for (double i = 3.0; i <= 9.0; i += 0.1)
// 	{
// 		fprintf_s(writer, "%lf,%lf\n", i, sp1.interpolate(i));
// 	}
// 	fclose(writer);
	quadratic_spline sp2(X, Y);
	result = sp2.interpolate(5.0);
	C = sp2.coeff();
	C.cout();
// 	fopen_s(&writer, "D:\\quadratic_spline.csv", "w");
// 	for (double i = 3.0; i <= 9.0; i += 0.1)
// 	{
// 		fprintf_s(writer, "%lf,%lf\n", i, sp2.interpolate(i));
// 	}
// 	fclose(writer);
	cubic_spline sp3(X, Y, notaknot());
	result = sp3.interpolate(5.0);
	C = sp3.coeff();
	C.cout();
// 	fopen_s(&writer, "D:\\cubic_spline.csv", "w");
// 	for (double i = 3.0; i <= 9.0; i += 0.1)
// 	{
// 		fprintf_s(writer, "%lf,%lf\n", i, sp3.interpolate(i));
// 	}
// 	fclose(writer);
	hermite hm(X, Y);
	result = hm.interpolate(5.0);
	C = hm.coeff();
	C.cout();
// 	fopen_s(&writer, "D:\\hermite_spline.csv", "w");
// 	for (double i = 3.0; i <= 9.0; i += 0.1)
// 	{
// 		fprintf_s(writer, "%lf,%lf\n", i, hm.interpolate(i));
// 	}
// 	fclose(writer);
	akima ak(X, Y);
	result = ak.interpolate(5.0);
	C = ak.coeff();
	C.cout();
// 	fopen_s(&writer, "D:\\akima_spline.csv", "w");
// 	for (double i = 3.0; i <= 9.0; i += 0.1)
// 	{
// 		fprintf_s(writer, "%lf,%lf\n", i, ak.interpolate(i));
// 	}
// 	fclose(writer);
}