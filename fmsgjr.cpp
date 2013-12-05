// fmsvalue.cpp - test fmsvalue
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include "fmsgjr.h"

using namespace fms;
static std::default_random_engine dre;

void test_bell()
{
	using polynomial::Bell;

	double x[] = {1,2,3,4};

	ensure (Bell(0,4,x) == 1);
	ensure (Bell(1,4,x) == x[0]);
	ensure (Bell(2,4,x) == x[0]*x[0] + x[1]);
	ensure (Bell(3,4,x) == x[0]*x[0]*x[0] + 3*x[0]*x[1] + x[2]);

	ensure (Bell(0,4,x,true) == 1);
	ensure (Bell(1,4,x,true) == x[0]);
	ensure (2*Bell(2,4,x,true) == x[0]*x[0] + x[1]);
	ensure (3*2*Bell(3,4,x,true) == x[0]*x[0]*x[0] + 3*x[0]*x[1] + x[2]);
}

void test_hermite()
{
}

void test_polynomial()
{
	test_bell();
	test_hermite();
}

void test_normal()
{
	double max, min, eps;
	eps = std::numeric_limits<double>::epsilon();

	max = -std::numeric_limits<double>::max();
	min = std::numeric_limits<double>::max();

	// 
	double ed_66[] = {
		-1.000000000000000, 
		-0.9999999999984625, 
		-0.99999998458274210, 
		-0.99997790950300141, 
		-0.99532226501895273, 
		-0.84270079294971487, 
		0, 0.84270079294971487, 
		0.99532226501895273, 
		0.99997790950300141, 
		0.99999998458274210, 
		0.9999999999984625, 
		1.000000000000000
	};
	for (int x = -6; x <= 6; x += 1) {
		double dy, y, y_;

		y = erf(x);
		y_ = ed_66[x + 6];
		dy = y - y_;
		max = (std::max<double>)(max, dy);
		min = (std::min<double>)(min, dy);
	}
	ensure (max < eps && min > -eps);

	// Table[N[Erfc[-x/Sqrt[2]]/2,17], {x, -6, 6}]
	// normaldist(x) x = -6,...,6
	double nd_66[] = {
		9.8658764503769814e-10,
		2.8665157187919391E-7,
		0.000031671241833119921,
		0.0013498980316300945,
		0.022750131948179207,
		0.15865525393145705,
		0.50000000000000000,
		0.84134474606854295,
		0.97724986805182079,
		0.99865010196836991,
		0.99996832875816688,
		0.99999971334842812,
		0.99999999901341235
	};

	max = -std::numeric_limits<double>::max();
	min = std::numeric_limits<double>::max();
	for (int x = -6; x <= 6; x += 1) {
		double dy, y, y_;

		y = distribution::standard_normal<>::cdf(x);
		y_ = nd_66[x + 6];
		dy = y - y_;
		max = (std::max<double>)(max, dy);
		min = (std::min<double>)(min, dy);
	}
	ensure (max < eps && min > -eps);
}
void test_distribution()
{
	test_normal();
}

template<class T>
void test_value()
{
	T f, s, k, t;

	f = 100;
	s = 0.2;
	k = 100;
	t = 0.25;

	model::black<> mb(f,s);

	T v, v_;
	v = value<T>(mb, instrument::option<>(k,t));
	v_ = value<T>(mb, instrument::put<>(k,t));
	ensure (v == v_);
	v = value<T>(mb, instrument::call<>(k,t));
	ensure (v == v_);
	v_ = value<T>(model::gjr<>(f,s), instrument::option<>(k,t));
	model::gjr<> mg(f,s);
	instrument::option<> o(k,t);
	v = value<T>(mg, o);

	model::black<float,int> mbfi(100., 1);
	v = value<T>(mb, instrument::option<>(k,t));
}

int main()
{
	try {
		test_polynomial();
		test_distribution();
		test_value<double>();
	}
	catch(const std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}

	return 0;
}