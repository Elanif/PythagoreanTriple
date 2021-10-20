#include<iostream>
#include<string>
#include<limits>
#include<cstdint>
#include<intrin.h>
#include<functional>

#include"SqrtApproximator.hpp"
#include"Constants.hpp"


void sum_and_subtraction_optimized(const bool trueint_falsedouble, int_fast64_t number) {
    if (trueint_falsedouble) {
        const int_fast64_t starting_number = 3;
        const int_fast64_t starting_sum = starting_number * starting_number;
        int_fast64_t hypothenuse_squared = number * number;
        int_fast64_t B = number - 1;
        int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum + B * B - number*number;
        /*uint_fast64_t max_iterations = 0;*/
        if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
        for (int_fast64_t A = (starting_number + 1) * 2 - 1; (A + 1) / 2 <= number / sqrt_2;) {
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
        std::cout << "sum_and_subtraction_optimized_int out\n\r";
    }
    else {
        const int_fast64_t starting_number = 3;
        const int_fast64_t starting_sum = starting_number * starting_number;
        int_fast64_t hypothenuse_squared = number * number;
        double B = number - 1;
        double sum = ((starting_sum + 1) - number) - number;//= starting_sum + B * B - number*number;
        /*uint_fast64_t max_iterations = 0;*/
        if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
        for (double A = (starting_number + 1) * 2 - 1; (A + 1) / 2 <= number / sqrt_2;) {
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
                    std::cout << "(" << int_fast64_t(A + 1) / 2 << "," << int_fast64_t(B) << "," << number << ") ";
                else
                    std::cout << "Error!\n\r";
            }

            A += 2;
        }
        std::cout << "sum_and_subtraction_optimized_double out\n\r";
    }
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
    multiplication_floor, //floor
    multiplication_ceil //ceil
};

int_fast64_t ceil_sqrt(int_fast64_t n) {
    return std::ceil(std::sqrt(n));
}

//sqrt_approx:N+->N+ such as 0<=sqrt_approx(x)<=ceil(sqrt(x))
void find_next_A(int_fast64_t& A, int_fast64_t& sum, const FindNextA& find_type, const std::function<int_fast64_t(int_fast64_t)>& sqrt_approx = ceil_sqrt) {
    //switch is so ugly in c++ and optimizer should do its thing
    if (find_type == FindNextA::square_root_1) {
        int_fast64_t temp = -A + sqrt_approx(A * A - sum); 
        if (temp <= 0) temp = 1;
        sum += temp * (2 * A + temp);
        A += temp;
    }
    else if (find_type == FindNextA::square_root_2) {
        int_fast64_t temp = sqrt_approx(A * A - sum);
        if (temp <= A) temp = A + 1;
        sum += temp * temp - A * A; //(temp-A)*(2*A-A+temp)
        A = temp;
    }
    else if (find_type == FindNextA::multiplication_floor || find_type == FindNextA::multiplication_ceil) {
        std::function<int_fast64_t(int_fast64_t, int_fast64_t)> division;
        if (find_type == FindNextA::multiplication_floor) division = [=](int_fast64_t a, int_fast64_t b)->int_fast64_t {
            return a / b;
        };
        if (find_type == FindNextA::multiplication_ceil) division = [=](int_fast64_t a, int_fast64_t b)->int_fast64_t {
            return a / b + ((a % b) != 0);
        };

        if (sum == 2 * A - 1) {
            sum = 0;
            ++A;
        }
        else {
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

void sum_and_subtraction_jump(int_fast64_t number, FindNextA method, const std::function<int_fast64_t(int_fast64_t)>& sqrt_approx = ceil_sqrt) {
    const int_fast64_t starting_number = 3;
    const int_fast64_t starting_sum = starting_number * starting_number;
    int_fast64_t hypothenuse_squared = number * number;
    int_fast64_t B = number - 1;
    int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum * starting_sum  + B * B - number*number;
    /*uint_fast64_t max_iterations = 0;*/
    if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
    for (int_fast64_t A = starting_number; A <= number / sqrt_2;) {

        find_next_A(A, sum, method, sqrt_approx);

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

    tester(std::cout, 0, 256);

    sum_and_subtraction(13);
    sum_and_subtraction_optimized(true, 13);
    sum_and_subtraction_jump(13, FindNextA::square_root_1);
    sum_and_subtraction_jump(13, FindNextA::multiplication_floor);
    sum_and_subtraction_jump(13, FindNextA::multiplication_ceil);

    int_fast64_t test_number = 5*(2LLu << 26 + 1);// 2LLu << 25 + 7;
    sum_and_subtraction(test_number);
    sum_and_subtraction_optimized(true, test_number);
    sum_and_subtraction_optimized(false, test_number);
    fast_ceil_sqrt_double sqrt_double;
    fast_ceil_sqrt_bitshift_lowacc sqrt_bitshift;
    sum_and_subtraction_jump(test_number, FindNextA::square_root_2);
    sum_and_subtraction_jump(test_number, FindNextA::square_root_2, [&](int_fast64_t n) {return sqrt_double(n); });
    sum_and_subtraction_jump(test_number, FindNextA::square_root_2, [&](int_fast64_t n) {return sqrt_bitshift(n); });
}
