#include<iostream>
#include<string>
#include<limits>
#include<cstdint>
const double sqrt_2 = std::sqrt(2.);

void sum_and_subtraction(int_fast64_t number) {
    const int_fast64_t starting_number = 3;
    const int_fast64_t starting_sum = starting_number * starting_number;
    int_fast64_t hypothenuse_squared = number * number;
    int_fast64_t B = number - 1;
    int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum + B * B - number*number;
    /*uint_fast64_t max_iterations = 0;*/
    if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
    for (int_fast64_t A = starting_number + 1; A <= number / sqrt_2; ++A) {
        sum += A * 2 - 1;

        /*uint_fast64_t iterations = 0;*/
        while (sum > 0 && B > 0) {
            /*++iterations;*/
            sum -= --B * 2 + 1;
        }
        /*if (iterations > max_iterations) max_iterations = iterations;*/

        if (sum == 0) {
            if (A * A + B * B == hypothenuse_squared)
                std::cout << "(" << A << "," << B << "," << number << ") ";
            else
                std::cout << "Error!\n\r";
        }
    }
    std::cout << "sum_and_subtraction out\n\r";
    /*std::cout << "max iterations = " << max_iterations << "\n\r";*/
}

enum class FindNextA {
    square_root_1,
    square_root_2,
    multiplication_1,
    multiplication_2,
};
void find_next_A(int_fast64_t& A, int_fast64_t& sum, const FindNextA& find_type) {
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
            int_fast64_t temp = (-sum + 2 * A - 2) / (2 * A - 1);
            if (temp == 0) temp = 1;
            sum += temp * (2 * A + temp);
            A += temp;
        }
    }
    else if (find_type == FindNextA::multiplication_2) {
        if (sum == 2 * A - 1) {
            sum = 0;
            ++A;
        }
        else {
            int_fast64_t temp = sum / (2 * A - 1);
            if (temp == 0) temp = 1;
            sum += temp * (2 * A + temp);
            A += temp;
        }
    }
}

void sum_and_subtraction_opt(int_fast64_t number) {
    const int_fast64_t starting_number = 3;
    const int_fast64_t starting_sum = starting_number * starting_number;
    int_fast64_t hypothenuse_squared = number * number;
    int_fast64_t B = number - 1;
    int_fast64_t sum = ((starting_sum + 1) - number) - number;//= starting_sum + B * B - number*number;
    /*uint_fast64_t max_iterations = 0;*/
    if (sum == 0) std::cout << "(" << starting_number << "," << B << "," << number << ") ";
    for (int_fast64_t A = starting_number; A <= number / sqrt_2;) {

        find_next_A(A, sum, FindNextA::multiplication_1);        

        /*uint_fast64_t iterations = 0;*/
        while (sum > 0 && B > 0) {
            /*++iterations;*/
            sum -= --B * 2 + 1;
        }
        /*if (iterations > max_iterations) max_iterations = iterations;*/

        if (sum == 0) {
            if (A * A + B * B == hypothenuse_squared)
                std::cout << "(" << A << "," << B << "," << number << ") ";
            else
                std::cout << "Error!\n\r";
        }
    }
    std::cout << "sum_and_subtraction_opt out\n\r";
    /*std::cout << "max iterations = " << max_iterations << "\n\r";*/
}

int main()
{
    enum class ProblemType {
        Cathetus,
        Hypothenuse,
        Undefined
    } problem_type = ProblemType::Undefined;

    std::string pb_type_string = "";
    while (problem_type == ProblemType::Undefined) {
        std::cout << "Cathetus or Hypothenuse? c/h ";
        std::cin >> pb_type_string;
        if (pb_type_string.length() > 0) {
            if (pb_type_string[0] == 'c' || pb_type_string[0] == 'C')
                problem_type = ProblemType::Cathetus;
            else if (pb_type_string[0] == 'h' || pb_type_string[0] == 'H')
                problem_type = ProblemType::Hypothenuse;
        }
    }

    if (problem_type == ProblemType::Hypothenuse) std::cout << "Hypothenuse\n\r";
    if (problem_type == ProblemType::Cathetus) std::cout << "Cathetus\n\r";

    std::cout << "Enter an integer ";
    int_fast64_t number = 0;
    std::cin >> number;
    std::cout << "trying with " << number << "\n\r";
    //if (number> std::numeric_limits<int_fast64_t>::min())
    number = 65;// (2llu << 30) - 5;
    std::cout << "trying with " << number << "\n\r";

    if (problem_type == ProblemType::Hypothenuse) {
        sum_and_subtraction(number);
        sum_and_subtraction_opt(number);
    }
    return 0;
}
