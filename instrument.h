// instrument.h
#pragma once

namespace fms {
	namespace instrument {

		template<class K = double, class T = double>
		struct option {
			K strike;
			T expiration;
			option(K k = 0, T t = 0)
				: strike(k), expiration(t)
			{ }
		};
		template<class K = double, class T = double>
		struct call : public option<K,T> {
			call(K k = 0, T t = 0)
				: option(k,t)
			{
				ensure (k >= 0);
			}
		};
		template<class K = double, class T = double>
		struct put : public option<K,T> {
			put(K k = 0, T t = 0)
				: option(k,t)
			{
				if (k >= 0)
					strike = -strike;
			}
		};

		// pw linear???

	} // instrument
} // fms