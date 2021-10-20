#pragma once
#include<intrin.h>
#include<cstdint>
#include<cmath>
#include"Constants.hpp"
#include<ostream>
#include<iostream>

const long double error_epsilon = 1.L - 1e-7L;

class fast_ceil_sqrt_bitshift_notwork { //doesnt always work with number with an odd number of bits (e.g. 1110000 = 120-> 1110 = 14) 
public:
    fast_ceil_sqrt_bitshift_notwork() {
    }
    int_fast64_t operator()(int_fast64_t n) {
        //if (n == 0) return 0;
        //if (n == 1) return 1;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        n = n >> (first_significant_bit >> 1);
        if (!(first_significant_bit & 1)) return n; //wrong one
        return n >> 1;
    }
};

class fast_ceil_sqrt_bitshift_lowacc { //faster? than fast_ceil_sqrt_double
public:
    fast_ceil_sqrt_bitshift_lowacc() {
    }
    int_fast64_t operator()(int_fast64_t n) {
        //if (n == 0) return 0;
        //if (n == 1) return 1;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        return n >> (1 + (first_significant_bit >> 1));
    }
};

class fast_ceil_sqrt_bitshift { //faster? than fast_ceil_sqrt_double, less accurate tho?
public:
    fast_ceil_sqrt_bitshift() {
    }
    int_fast64_t operator()(int_fast64_t n) {
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        n = n >> (1+(first_significant_bit >> 1));
        if (first_significant_bit & 1) return n;
        else {
            const int_fast64_t mask = 0b1010101010101010101010101010101010101010101010101010101010101111ll;
            return n & mask;
        }
    }
};

//grabs the first 2 significant bits in a number e.g.: 101110 -> 101000 = 2^6+2^4 = 2^2(1+2^4) -> sqrt(101110) ~ sqrt(2^2)*sqrt(1+2^4)
class fast_ceil_sqrt_double { //TODO what if ceil fails? floating point fails anyway at big numbers, so at some point it would be wise to multiply everything by like (1-10^-7)
public:
    long double first_exponent_root[64];
    long double second_exponent_root[64];
    int_fast64_t result[63][63];
    fast_ceil_sqrt_double() {
        first_exponent_root[0] = 1; //careful, it's not shifted
        second_exponent_root[0] = 1; //when grabbing the significant bit if a number is 0 or 1 the result is still 0, so if the bumber is 1 we instead count it as 1, the array is shifted
        for (int_fast64_t i = 1; i < 64; ++i) {
            first_exponent_root[i] = std::sqrtl(std::powl(2, i));
            second_exponent_root[i] = std::sqrtl(1 + std::powl(2, i - 1));
        }
        for (int_fast64_t i = 0; i < 63; ++i) {
            for (int_fast64_t j = 0; j < 63; ++j) {
                result[i][j] = 0;
            }
        }
        for (int_fast64_t i = 0; i < 63; ++i) {
            for (int_fast64_t j = 0; j <= i; ++j) {
                if (j == 0) {
                    //all of the next in this if can be handled by result[i][j]=std::ceill(first_exponent_root[i]);
                    if (i == 0) result[i][j] = 1;
                    else if (i % 2 == 0)
                        result[i][j] = int_fast64_t(2) << (i / 2 -1);
                    else if (i>=2)
                        result[i][j] = std::ceil(error_epsilon*((int_fast64_t(2) << (i / 2 -1)) * sqrtl(2)));
                    else
                        result[i][j] = std::ceil(sqrtl(2));
                }
                else {
                    if (j == 1) //not necessary
                        result[i][j] = std::ceill(error_epsilon*second_exponent_root[i + 1]);
                    else
                        result[i][j] = std::ceill(error_epsilon*first_exponent_root[j - 1] * second_exponent_root[i - j + 2]);

                    //failsafe
                    for (int_fast64_t m = 0; m <= i - 1; ++m)
                        if (result[i - 1][m] > result[i][j]) result[i][j] = result[i - 1][m];
                }
            }
        }
    }
    int_fast64_t operator()(int_fast64_t n) {
        if (n == 0) return 0;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        unsigned long second_significant_bit = 0;
        int_fast64_t removed_bit = n & (~(int_fast64_t(1) << first_significant_bit));
        //int_fast64_t removed_bit = n - (int_fast64_t(1) << first_significant_bit);
        if (removed_bit) {
            _BitScanReverse64(&second_significant_bit, removed_bit);
            ++second_significant_bit;
        }
        return result[first_significant_bit][second_significant_bit];
    }
};


class fast_ceil_sqrt_bitshift_ceil { //doesnt work for obvious reasons //also obsoleted by fast_ceil_sqrt_double
public:
    fast_ceil_sqrt_bitshift_ceil() {
    }
    int_fast64_t operator()(int_fast64_t n) {
        if (n == 0) return 0;
        if (n == 1) return 1;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        unsigned long least_significant_bit = 0;
        _BitScanForward64(&least_significant_bit, n);
        first_significant_bit = first_significant_bit / 2 + first_significant_bit & 1;
        return n >> first_significant_bit;
    }
};



class fast_ceil_sqrt_int { //obsoleted by fast_ceil_sqrt_double
public:
    int_fast64_t first_exponent_root[64];
    int_fast64_t second_exponent_root[64];
    fast_ceil_sqrt_int() {
        first_exponent_root[0] = 0;
        second_exponent_root[0] = 1;
        for (int_fast64_t i = 1; i < 64; ++i) {
            first_exponent_root[i] = std::sqrt(std::pow(2, i - 1));
            second_exponent_root[i] = std::sqrt(1 + std::pow(2, i - 1));
        }
    }
    int_fast64_t operator()(int_fast64_t n) {
        if (n == 0) return 0;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        unsigned long second_significant_bit = 0;
        int_fast64_t removed_bit = n & (~(int_fast64_t(1) << (first_significant_bit)));
        if (removed_bit != 0) {
            _BitScanReverse64(&second_significant_bit, removed_bit);
        }
        else return first_exponent_root[first_significant_bit + 1];

        return first_exponent_root[second_significant_bit + 1] * second_exponent_root[first_significant_bit - second_significant_bit + 1];
    }
};

class fast_ceil_sqrt_low_high { //doesnt work for obvious reasons //also obsoleted by fast_ceil_sqrt_double
public:
    int_fast64_t first_exponent_root_low[64];
    int_fast64_t first_exponent_root_high[64];

    int_fast64_t second_exponent_root_low[64];
    int_fast64_t second_exponent_root_high[64];
    fast_ceil_sqrt_low_high() {
        /*first_exponent_root_low[0] = 0;
        first_exponent_root_low[1] = 1;
        first_exponent_root_high[0] = 0;
        first_exponent_root_high[1] = 1;
        for (int_fast64_t i = 3; i < 64; i += 2)
            first_exponent_root_low[i] = first_exponent_root_high[i] = first_exponent_root_low[i - 2] << 1;
        for (int_fast64_t i = 2; i < 64; i += 2) {
            first_exponent_root_low[i] = std::sqrt(int_fast64_t(2) << (i - 1));
            first_exponent_root_high[i] = std::ceil(std::sqrt(int_fast64_t(2) << (i - 1)));
        }
        second_exponent_root_low[0] = 1;
        second_exponent_root_low[1] = 1;
        second_exponent_root_high[0] = 1;
        second_exponent_root_high[1] = 2;
        for (int_fast64_t i = 2; i < 64; ++i) {
            second_exponent_root_low[i] = std::sqrt(1 + (int_fast64_t(2) << (i - 1)));
            second_exponent_root_high[i] = std::ceil(std::sqrt(1 + (int_fast64_t(2) << (i - 1))));
        }*/
        first_exponent_root_low[0] = 1;
        first_exponent_root_high[0] = 1;
        second_exponent_root_low[0] = 1;
        second_exponent_root_high[0] = 1;
        for (int_fast64_t i = 1; i < 64; ++i) {
            first_exponent_root_low[i] = std::sqrt(std::pow(2, i - 1));
            first_exponent_root_high[i] = std::ceil(std::sqrt(std::pow(2, i - 1)));
            second_exponent_root_low[i] = std::sqrt(1 + std::pow(2, i - 1));
            second_exponent_root_high[i] = std::ceil(std::sqrt(1 + std::pow(2, i - 1)));
        }
    }
    int_fast64_t operator()(int_fast64_t n) {
        if (n == 0) return 0;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        unsigned long second_significant_bit = 0;
        int_fast64_t removed_bit = n & (~(int_fast64_t(1) << (first_significant_bit)));
        if (removed_bit != 0) {
            _BitScanReverse64(&second_significant_bit, removed_bit);
        }
        else return first_exponent_root_high[first_significant_bit + 1];
        if ((second_significant_bit % 2) == 1) {
            return first_exponent_root_high[second_significant_bit + 1] * second_exponent_root_low[first_significant_bit - second_significant_bit + 1];
        }
        else {
            return first_exponent_root_low[second_significant_bit + 1] * second_exponent_root_high[first_significant_bit - second_significant_bit + 1];
        }
    }
};

class fast_ceil_sqrt_high { //doesnt work for obvious reasons //also obsoleted by fast_ceil_sqrt_double
public:
    int_fast64_t first_exponent_root_low[64];
    int_fast64_t first_exponent_root_high[64];

    int_fast64_t second_exponent_root_low[64];
    int_fast64_t second_exponent_root_high[64];
    fast_ceil_sqrt_high() {
        first_exponent_root_low[0] = 1;
        first_exponent_root_high[0] = 1;
        second_exponent_root_low[0] = 1;
        second_exponent_root_high[0] = 1;
        for (int_fast64_t i = 1; i < 64; ++i) {
            first_exponent_root_low[i] = std::sqrt(std::pow(2, i - 1));
            first_exponent_root_high[i] = std::ceil(std::sqrt(std::pow(2, i - 1)));
            second_exponent_root_low[i] = std::sqrt(1 + std::pow(2, i - 1));
            second_exponent_root_high[i] = std::ceil(std::sqrt(1 + std::pow(2, i - 1)));
        }
    }
    int_fast64_t operator()(int_fast64_t n) {
        if (n == 0) return 0;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        unsigned long second_significant_bit = 0;
        int_fast64_t removed_bit = n & (~(int_fast64_t(1) << (first_significant_bit)));
        if (removed_bit != 0) {
            _BitScanReverse64(&second_significant_bit, removed_bit);
        }
        else return first_exponent_root_high[first_significant_bit + 1];
        return std::max(first_exponent_root_high[second_significant_bit + 1] * second_exponent_root_low[first_significant_bit - second_significant_bit + 1], first_exponent_root_low[second_significant_bit + 1] * second_exponent_root_high[first_significant_bit - second_significant_bit + 1]);
    }
};

class fast_ceil_sqrt_low { //doesnt work for obvious reasons //also obsoleted by fast_ceil_sqrt_double
public:
    int_fast64_t first_exponent_root_low[64];
    int_fast64_t first_exponent_root_high[64];

    int_fast64_t second_exponent_root_low[64];
    int_fast64_t second_exponent_root_high[64];
    fast_ceil_sqrt_low() {
        first_exponent_root_low[0] = 1;
        first_exponent_root_high[0] = 1;
        second_exponent_root_low[0] = 1;
        second_exponent_root_high[0] = 1;
        for (int_fast64_t i = 1; i < 64; ++i) {
            first_exponent_root_low[i] = std::sqrt(std::pow(2, i - 1));
            first_exponent_root_high[i] = std::ceil(std::sqrt(std::pow(2, i - 1)));
            second_exponent_root_low[i] = std::sqrt(1 + std::pow(2, i - 1));
            second_exponent_root_high[i] = std::ceil(std::sqrt(1 + std::pow(2, i - 1)));
        }
    }
    int_fast64_t operator()(int_fast64_t n) {
        if (n == 0) return 0;
        unsigned long first_significant_bit = 0;
        _BitScanReverse64(&first_significant_bit, n);
        unsigned long second_significant_bit = 0;
        int_fast64_t removed_bit = n & (~(int_fast64_t(1) << (first_significant_bit)));
        if (removed_bit != 0) {
            _BitScanReverse64(&second_significant_bit, removed_bit);
        }
        else return first_exponent_root_high[first_significant_bit + 1];
        return std::min(first_exponent_root_high[second_significant_bit + 1] * second_exponent_root_low[first_significant_bit - second_significant_bit + 1], first_exponent_root_low[second_significant_bit + 1] * second_exponent_root_high[first_significant_bit - second_significant_bit + 1]);
    }
};

template<class T>
void tester(std::basic_ostream<T>& outstream, int_fast64_t min, int_fast64_t max) {


    /*fast_ceil_sqrt_low_high sqrt_low_high;
    fast_ceil_sqrt_low sqrt_low;
    fast_ceil_sqrt_high sqrt_high;
    fast_ceil_sqrt_int sqrt_int;*/
    fast_ceil_sqrt_double sqrt_double;
    //fast_ceil_sqrt_bitshift_notwork sqrt_bitshift_notwork;
    fast_ceil_sqrt_bitshift sqrt_bitshift;
    fast_ceil_sqrt_bitshift_lowacc sqrt_bitshift_lowacc;

    sqrt_double(std::numeric_limits<int_fast64_t>::max());

    for (int_fast64_t i = min; i < max; ++i) {
        int_fast64_t real = std::ceil(std::sqrt(i));
        /*int_fast64_t low_high = sqrt_low_high(i);
        int_fast64_t low = sqrt_low(i);
        int_fast64_t high = sqrt_high(i);
        int_fast64_t int_ = sqrt_int(i);*/
        int_fast64_t double_ = sqrt_double(i);
        //int_fast64_t bitshift_notwork_ = sqrt_bitshift_notwork(i);
        int_fast64_t bitshift_ = sqrt_bitshift(i);
        int_fast64_t bitshift_lowacc_ = sqrt_bitshift_lowacc(i);
        outstream << "i = " << i << " real= " << real << " ";
        /*outstream << "low_high result = " << low_high << " ";
        outstream << "low result = " << low << " ";
        outstream << "high result = " << high << " ";
        outstream << "int_ result = " << int_ << " ";*/
        outstream << "double_= " << double_ << " ";
        outstream << "bitshift_= " << bitshift_ << " ";
        outstream << "lowacc_= " << bitshift_lowacc_ << "\n\r";
    }
}