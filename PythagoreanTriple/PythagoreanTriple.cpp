#include"Constants.hpp"

#include"IterativePythagoreanFinders.hpp"
#include"SqrtApproximator.hpp"
#include"FactoringAlgorithms.hpp"


int main()
{
    std::cout.sync_with_stdio(false);
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

    if (problem_type == ProblemType::Hypothenuse) {
        sum_and_subtraction(number);
        sum_and_subtraction_optimized(true, number);
        fast_ceil_sqrt_double sqrt_double;
        //fast_ceil_sqrt_bitshift_lowacc sqrt_bitshift;
        sum_and_subtraction_jump(number, FindNextA::square_root_2);
    }

    /*tester(std::cout, 0, 256);

    sum_and_subtraction(13);
    sum_and_subtraction_optimized(true, 13);
    sum_and_subtraction_jump(13, FindNextA::square_root_1);*/
    //sum_and_subtraction_jump(13, FindNextA::multiplication_floor);
    //sum_and_subtraction_jump(13, FindNextA::multiplication_ceil);

    //int_fast64_t test_number = 5ll * 13ll * 17ll * 29ll * 37ll * 43ll;//5*(2LLu << 26 + 1);
    //sum_and_subtraction(test_number);
    //sum_and_subtraction_optimized(true, test_number);
    //sum_and_subtraction_optimized(false, test_number);
    //fast_ceil_sqrt_double sqrt_double;
    //fast_ceil_sqrt_bitshift_lowacc sqrt_bitshift;
    //sum_and_subtraction_jump(test_number, FindNextA::square_root_2);
    //sum_and_subtraction_jump(test_number, FindNextA::square_root_2, [&](int_fast64_t n) {return sqrt_double(n); });
    //sum_and_subtraction_jump(test_number, FindNextA::square_root_2, [&](int_fast64_t n) {return sqrt_bitshift(n); });
    std::cout << brent(1371839117) << "\n\r";
    std::cout << std::gcd(3824685, 5422289) << "\n\r";
    factor(1371839117);
}
