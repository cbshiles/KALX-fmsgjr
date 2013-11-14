// fmsvalue.cpp - test fmsvalue
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
//#include "exp.h"
//#include "fms_distribution.h"
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
		test_value<double>();
	}
	catch(const std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}

	return 0;
}