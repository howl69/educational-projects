#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define X1_A -5
#define X1_B 5
#define X2_A -10
#define X2_B -0.01
#define X3_A -10
#define X3_B -0.01
#define EPS_INTEGRAL 0.001/4
#define EPS_ROOT 0.000006
#define TEST_EPS_INTEGRAL 0.00001

extern long long int iterations;

extern double root(double func1(double), double func2(double), double a, double b, double eps);
extern double integral(double func(double), double a, double b, double n);
extern double test_root_eps();
extern void test_root();
extern void test_integral();
extern double f1(double x);
extern double f2(double x);
extern double f3(double x);
extern double df1(double x);
extern double df2(double x);
extern double df3(double x);
extern double f4(double x);
extern double f5(double x);
extern double f6(double x);
extern double f7(double x);
extern double f8(double x);
extern double f9(double x);

int main(int argc, char **argv)
{
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "-help") == 0)
        {
            printf("-help - all keys\n");
            printf("-roots - roots\n");
            printf("-iterations - all iterations of the equation solution\n");
            printf("-root_test - test root of function\n");
            printf("-integral_test - test integral of function\n");
            printf("-example_root_test - root test of functions:\n sin(x)=cos(x)\n 2^x+2=4^x\n log2(x-1)=-log2(x+1)\n");
            printf("-example_integral_test - integral test of functions:\n sin(x), on [0, pi]\n 2^x+2, on [-10, 10]\n 5/x^2, on [1.1, 10]\n");
            break;
        }
    double test;
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "-example_root_test") == 0)
        {
            test_root();
            break;
        }
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "-example_integral_test") == 0)
        {
            test_integral();
            break;
        }
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-root_test") == 0)
        {
            double (*func1)(double), (*func2)(double);
            if (atoi(argv[i + 1]) == 1)
                func1 = f1;
            if (atoi(argv[i + 1]) == 2)
                func1 = f2;
            if (atoi(argv[i + 1]) == 3)
                func1 = f3;
            if (atoi(argv[i + 1]) == 4)
                func1 = f4;
            if (atoi(argv[i + 1]) == 5)
                func1 = f5;
            if (atoi(argv[i + 1]) == 6)
                func1 = f6;
            if (atoi(argv[i + 1]) == 7)
                func1 = f7;
            if (atoi(argv[i + 1]) == 8)
                func1 = f8;
            if (atoi(argv[i + 1]) == 9)
                func1 = f9;
            if (atoi(argv[i + 2]) == 1)
                func2 = f1;
            if (atoi(argv[i + 2]) == 2)
                func2 = f2;
            if (atoi(argv[i + 2]) == 3)
                func2 = f3;
            if (atoi(argv[i + 2]) == 4)
                func2 = f4;
            if (atoi(argv[i + 2]) == 5)
                func2 = f5;
            if (atoi(argv[i + 2]) == 6)
                func2 = f6;
            if (atoi(argv[i + 2]) == 7)
                func2 = f7;
            if (atoi(argv[i + 2]) == 8)
                func2 = f8;
            if (atoi(argv[i + 2]) == 9)
                func2 = f9;
            if (strcmp(argv[i + 1], "d1") == 0)
                func1 = df1;
            if (strcmp(argv[i + 2], "d1") == 0)
                func2 = df1;
            if (strcmp(argv[i + 1], "d2") == 0)
                func1 = df2;
            if (strcmp(argv[i + 2], "d2") == 0)
                func2 = df2;
            if (strcmp(argv[i + 1], "d3") == 0)
                func1 = df3;
            if (strcmp(argv[i + 2], "d3") == 0)
                func2 = df3;
            test = root(func1, func2, atof(argv[i + 3]), atof(argv[i + 4]), atof(argv[i + 5]));
        printf("Root test result = %.10lf\n", test);
        break;
        }
    }
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-integral_test") == 0)
        {
            if (atoi(argv[i + 2]) == 1)
                test = integral(f1, atof(argv[i + 3]), atof(argv[i + 4]), atoi(argv[i + 5]));
            if (atoi(argv[i + 1]) == 2)
                test = integral(f2, atof(argv[i + 3]), atof(argv[i + 4]), atoi(argv[i + 5]));
            if (atoi(argv[i + 1]) == 3)
                test = integral(f3, atof(argv[i + 3]), atof(argv[i + 4]), atoi(argv[i + 5]));
            printf("Integral test result= %.10lf\n", test);
            break;
        }
    }
    double x1 = root(f1, f2, X1_A, X1_B, EPS_ROOT);
    double x2 = root(f1, f3, X2_A, X2_B, EPS_ROOT);
    double x3 = root(f2, f3, X3_A, X3_B, EPS_ROOT);
    double s = 0;
    s += integral(f3, x2, x3, EPS_INTEGRAL);
    s += integral(f2, x3, x1, EPS_INTEGRAL);
    s -= integral(f1, x2, x1, EPS_INTEGRAL);
    printf("Results:\n");
    printf("S = %.5lf\n", s);
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "-roots") == 0)
        {
            printf("x1 = %.5lf\n", x1);
            printf("x2 = %.5lf\n", x2);
            printf("x3 = %.5lf\n", x3);
            break;
        }
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "-iterations") == 0)
        {
            printf("Iterations = %lld\n", iterations);
            break;
        }
    return 0;
}
