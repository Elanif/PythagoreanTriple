#pragma once
#include "Constants.hpp"
#include"random.hpp"

using Random = effolkronium::random_static;

typedef std::tuple<int_fast64_t, int_fast64_t> pair;
typedef std::tuple<int_fast64_t, int_fast64_t, int_fast64_t> triple; 

long long trial_division3(long long n, std::map< int_fast64_t, int_fast64_t>& previous, long long max) {
	auto add = [&](int_fast64_t d) {
		auto it = previous.find(d);
		if (it != previous.end())
			--it->second;
		else
			previous[d] = 1;
	};
	for (int d : {2, 3, 5}) {
		while (n % d == 0) {
			add(d);
			n /= d;
		}
	}
	static std::array<int, 8> increments = { 4, 2, 4, 2, 4, 6, 2, 6 };
	int i = 0;
	for (long long d = 7; d * d <= max; d += increments[i++]) {
		while (n % d == 0) {
			add(d);
			n /= d;
		}
		if (i == 8)
			i = 0;
	}
	return n;
}

void factor_recursive(int_fast64_t n, std::map< int_fast64_t, int_fast64_t>& previous, int_fast64_t depth) {
	if (n == 1) return;
	auto x0 = Random::get();
	auto c = Random::get();
	int_fast64_t result = brent(n,x0,c);
	auto it = previous.find(result);
	if (it != previous.end()) //if it's in then either it's fucked or it's ok to do this
		it->second++;
	else {
		if (result == n || result == 1) {
			if (depth <= 0) {// give up
				previous[n]=1;
			}
			else { //try again
				factor_recursive(n, previous, depth - 1);
			}
		}
		else {
			factor_recursive(n / result, previous, depth + 1);
			factor_recursive(result, previous, depth + 1);
		}
	}

	/*if (depth <= 0) {
		if (result == n) {
			auto it = previous.find(result);
			if (it != previous.end())
				it->second++;
			else previous[result] = 1;
		}
		else {
			factor_recursive(n/result, previous, depth + 1);
		}
	}
	else {
		auto it = previous.find(result);
		if (it != previous.end()) {
			it->second++;
			if (result != n && result!=1) 
				factor_recursive(n / result, previous, depth + 1);
		}
		else {
			if (result!=n) factor_recursive(n / result, previous, depth + 1);
			if (result!=1) factor_recursive(result, previous, depth + 1);
		}
	}*/
	/*if (result == n) {
		if (depth <= 0) {
			auto it = previous.find(result);
			if (it != previous.end())
				it->second++;
			else previous[result] = 1;
		}
		else {
			factor_recursive(n, previous, depth - 1);
		}
	}
	else {
		auto it = previous.find(result);
		if (it != previous.end()) {
			it->second++;
			factor_recursive(n / result, previous, depth + 1);
		}
		else {
			previous[result] = 1;
			factor_recursive(n / result, previous, depth + 1);
			factor_recursive(result, previous, depth + 1);
		}
		
	}*/
}

std::vector<pair> factor(int_fast64_t n) {
	std::vector<pair> result;
	std::map< int_fast64_t, int_fast64_t> factorization_map;
	//factorization_map[1] = 1;
	n=trial_division3(n, factorization_map, 100000);
	factor_recursive(n, factorization_map, 10);
	return result;
};

std::tuple<pair> euclidian_algorithm(int_fast64_t big, int_fast64_t small, int_fast64_t cutoff) {
	return pair(1, 0);
};

std::vector<triple> generate_triples(const std::vector<pair>& factorization) {
	int_fast64_t multiplier = 1;
	std::vector<pair> good_primes;
	for (auto const& it : factorization) {
		if (std::get<0>(it) % 4 == 1) {
			good_primes.push_back(it);
		}
		else {
			multiplier *= ipow(std::get<1>(it),std::get<0>(it));
		}
	}
	std::vector<triple> result;
	return result;	
}

std::vector<triple> generate_triples(int_fast64_t n) {
	return generate_triples(factor(n));
}