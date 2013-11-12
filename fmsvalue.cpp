// fmsvalue.cpp - test fmsvalue
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
//#include "exp.h"
//#include "fms_distribution.h"
#include "fmsvalue.h"

using namespace fms;
static std::default_random_engine dre;

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
//	v_ = value<T>(model::jarrow_rudd<>(f,s), instrument::option<>(k,t));

	model::black<float,int> mbfi(100., 1);
	v = value<T>(mb, instrument::option<>(k,t));
}

int main()
{
	try {
		test_value<double>();
	}
	catch(const std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}

	return 0;
}