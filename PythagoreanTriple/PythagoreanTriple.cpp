#include<iostream>
#include<string>

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
    long number = 0;
    std::cin >> number;
    const long starting_number = 3;
    const long starting_sum = starting_number * starting_number;

    if (problem_type == ProblemType::Hypothenuse) {
        long hypothenuse_squared = number * number;
        long A = starting_number;
        long B = number - 1;
        long sum = starting_sum+B*B;
        if (sum == hypothenuse_squared) std::cout << "(" << A << "," << B << "," << number << ") ";
        for (long add = 1; add < number / 2; ++add) {
            //if ()
        }
    }

    return 0;
}
