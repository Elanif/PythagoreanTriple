#include<iostream>
#include<string>
#include<limits>
#include<cstdint>
#include<intrin.h>
#include"SqrtApproximator.hpp"
#include"Constants.hpp"

void sum_and_subtraction_optimized(int_fast64_t number) {
    const int_fast64_t starting_number = 3;
    const int_fast64_t starting_sum = starting_number * starting_number;
    int_fast64_t hypothenuse_squared = number * number;
    int_fast64_t B = number - 1;
    int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum + B * B - number*number;
    /*uint_fast64_t max_iterations = 0;*/
    if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
    for (int_fast64_t A = (starting_number + 1)*2-1; (A+1)/2 <= number / sqrt_2;) {
        sum += A;
        /*uint_fast64_t iterations = 0;*/
        if (sum > 0 && B > 0) {
            /*++iterations;*/
            sum -= --B * 2 + 1;
        }
        /*if (iterations > max_iterations) max_iterations = iterations;*/

        if (sum == 0) {
            uint_fast64_t a_ = (A + 1) / 2;
            uint_fast64_t b_ = B;
            uint_fast64_t c_ = number;
            if (a_ * a_ + b_ * b_ == c_ * c_) 
                std::cout << "(" << (A + 1) / 2 << "," << B << "," << number << ") ";
            else
                std::cout << "Error!\n\r";
        }

        A += 2;
    }
    std::cout << "sum_and_subtraction_optimized out\n\r";
    /*std::cout << "max iterations = " << max_iterations << "\n\r";*/
}

void sum_and_subtraction(int_fast64_t number) {
    const int_fast64_t starting_number = 3;
    const int_fast64_t starting_sum = starting_number * starting_number;
    int_fast64_t hypothenuse_squared = number * number;
    int_fast64_t B = number - 1;
    int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum + B * B - number*number;
    /*uint_fast64_t max_iterations = 0;*/
    if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
    for (int_fast64_t A = starting_number + 1; A <= number / sqrt_2;) {
        sum += A * 2 - 1;

        /*uint_fast64_t iterations = 0;*/
        if (sum > 0 && B > 0) {
            /*++iterations;*/
            sum -= --B * 2 + 1;
        }
        /*if (iterations > max_iterations) max_iterations = iterations;*/

        if (sum == 0) {
            uint_fast64_t a_ = A;
            uint_fast64_t b_ = B;
            uint_fast64_t c_ = number;
            if (a_ * a_ + b_ * b_ == c_ * c_) {
                std::cout << "(" << A << "," << B << "," << number << ") ";
            }
            else
                std::cout << "Error!\n\r";
        }

        ++A;
    }
    std::cout << "sum_and_subtraction out\n\r";
    /*std::cout << "max iterations = " << max_iterations << "\n\r";*/
}

enum class FindNextA {
    square_root_1,
    square_root_2,
    multiplication_1
};
enum class Approx {
    floor,
    ceil
};
void find_next_A(int_fast64_t& A, int_fast64_t& sum, const FindNextA& find_type, const Approx& approx) {
    //switch is so ugly in c++ and optimizer should do its thing
    if (find_type == FindNextA::square_root_1) {
        int_fast64_t temp = std::ceil(-A + std::sqrt(A * A - sum));
        if (temp == 0) temp = 1;
        sum += temp * (2 * A + temp);
        A += temp;
    }
    else if (find_type == FindNextA::square_root_2) {
        int_fast64_t temp = std::ceil(std::sqrt(A * A - sum));
        if (temp == A) temp = A + 1;
        sum += temp * temp - A * A; //(temp-A)*(2*A-A+temp)
        A = temp;
    }
    else if (find_type == FindNextA::multiplication_1) {
        if (sum == 2 * A - 1) {
            sum = 0;
            ++A;
        }
        else {
            auto division=[=](int_fast64_t a, int_fast64_t b)->int_fast64_t{
                if (approx == Approx::floor) return a / b;
                else if (approx == Approx::ceil) return a / b + ((a % b)!=0);
            };

            int_fast64_t temp = division(-sum, 2 * A - 1);
            if (temp == 0) temp = 1;
            else {
                int_fast64_t y = 2 * (temp + A) - 1;
                int_fast64_t temp2 = division(-sum, y); //approx (x*sqrt(2)+y)/(1+sqrt(2)) -> (x*19+y*13)/32
                temp = temp2;//temp = (temp * 19 + temp2 * 13) / 32;
            }

            sum += temp * (2 * A + temp);
            A += temp;
        }
    }
}

void sum_and_subtraction_jump(int_fast64_t number, FindNextA method, Approx approx) {
    const int_fast64_t starting_number = 3;
    const int_fast64_t starting_sum = starting_number * starting_number;
    int_fast64_t hypothenuse_squared = number * number;
    int_fast64_t B = number - 1;
    int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum * starting_sum  + B * B - number*number;
    /*uint_fast64_t max_iterations = 0;*/
    if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
    for (int_fast64_t A = starting_number; A <= number / sqrt_2;) {

        find_next_A(A, sum, method, approx);

        /*uint_fast64_t iterations = 0;*/
        if (sum > 0 && B > 0) {
            /*++iterations;*/
            sum -= --B * 2 + 1;
        }
        /*if (iterations > max_iterations) max_iterations = iterations;*/

        if (sum == 0 && A <= number / sqrt_2) {
            uint_fast64_t a_ = A;
            uint_fast64_t b_ = B;
            uint_fast64_t c_ = number;
            if (a_ * a_ + b_ * b_ == c_*c_) { 
                std::cout << "(" << A << "," << B << "," << number << ") ";
            }
            else
                std::cout << "Error!\n\r";
            sum -= --B * 2 + 1;
            sum += 2 * A + 1;
            ++A;

        }

    }
    std::cout << "sum_and_subtraction_jump out\n\r";
    /*std::cout << "max iterations = " << max_iterations << "\n\r";*/
}

int main()
{
    enum class ProblemType {
        Cathetus,
        Hypothenuse,
        Undefined
    } problem_type = ProblemType::Undefined;

    /*fast_ceil_sqrt_low_high sqrt_low_high;
    fast_ceil_sqrt_low sqrt_low;
    fast_ceil_sqrt_high sqrt_high;
    fast_ceil_sqrt_int sqrt_int;*/
    fast_ceil_sqrt_double sqrt_double;
    fast_ceil_sqrt_bitshift_notwork sqrt_bitshift_notwork;
    
    sqrt_double(std::numeric_limits<int_fast64_t>::max());

    for (int_fast64_t i = 0; i < 260; ++i) {
        int_fast64_t real= std::ceil(std::sqrt(i));
        /*int_fast64_t low_high = sqrt_low_high(i);
        int_fast64_t low = sqrt_low(i);
        int_fast64_t high = sqrt_high(i);
        int_fast64_t int_ = sqrt_int(i);*/
        int_fast64_t double_ = sqrt_double(i);
        int_fast64_t bitshift_notwork_ = sqrt_bitshift_notwork(i);
        std::cout << "Sqrt of " << i << " real result = " << real << " ";
        /*std::cout << "low_high result = " << low_high << " ";
        std::cout << "low result = " << low << " ";
        std::cout << "high result = " << high << " ";
        std::cout << "int_ result = " << int_ << " ";*/
        std::cout << "double_ result = " << double_ << " ";
        std::cout << "bitshift_ result = " << bitshift_notwork_ << "\n\r";
    }

    //std::string pb_type_string = "";
    //while (problem_type == ProblemType::Undefined) {
    //    std::cout << "Cathetus or Hypothenuse? c/h ";
    //    std::cin >> pb_type_string;
    //    if (pb_type_string.length() > 0) {
    //        if (pb_type_string[0] == 'c' || pb_type_string[0] == 'C')
    //            problem_type = ProblemType::Cathetus;
    //        else if (pb_type_string[0] == 'h' || pb_type_string[0] == 'H')
    //            problem_type = ProblemType::Hypothenuse;
    //    }
    //}

    //if (problem_type == ProblemType::Hypothenuse) std::cout << "Hypothenuse\n\r";
    //if (problem_type == ProblemType::Cathetus) std::cout << "Cathetus\n\r";

    //std::cout << "Enter an integer ";
    //int_fast64_t number = 0;
    //std::cin >> number;
    //std::cout << "trying with " << number << "\n\r";
    ////if (number> std::numeric_limits<int_fast64_t>::min())
    //number = 66000001;// (2llu << 30) - 5;
    //std::cout << "trying with " << number << "\n\r";

    //if (problem_type == ProblemType::Hypothenuse) {
    //    sum_and_subtraction(number);
    //    sum_and_subtraction_optimized(number);
    //    sum_and_subtraction_jump(number, FindNextA::multiplication_1, Approx::floor);
    //}

    sum_and_subtraction(13);
    sum_and_subtraction_optimized(13);
    sum_and_subtraction_jump(13, FindNextA::square_root_1, Approx::ceil);

    sum_and_subtraction(26111111);
    sum_and_subtraction_optimized(26111111);
    sum_and_subtraction_jump(26111111, FindNextA::multiplication_1, Approx::ceil);
}
