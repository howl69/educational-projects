#include<stdio.h>
#include<math.h>

#define X1_MIN 1
#define X1_MAX 2
#define X2_MIN -3
#define X2_MAX -2
#define X3_MIN -1
#define X3_MAX -0.01
#define N_START 10
#define CONST_P 15
#define X1_A -5
#define X1_B 5
#define X2_A -10
#define X2_B -0.01
#define X3_A -10
#define X3_B -0.01
#define TEST_ROOT_X1_A 0
#define TEST_ROOT_X1_B 3.1415926
#define TEST_ROOT_X2_A -10
#define TEST_ROOT_X2_B 10
#define TEST_ROOT_X3_A 1.1
#define TEST_ROOT_X3_B 10
#define TEST_INTEGRAL_S1_A 0
#define TEST_INTEGRAL_S1_B 3.1415926
#define TEST_INTEGRAL_S2_A 0
#define TEST_INTEGRAL_S2_B 2
#define TEST_INTEGRAL_S3_A 2
#define TEST_INTEGRAL_S3_B 4
#define EPS_ROOT_START 0.1
#define EPS_INTEGRAL 0.001/4
#define EPS_ROOT 0.000006
#define TEST_EPS_INTEGRAL 0.00001
#define EPS 0.0001

long long int iterations = 0;
extern double f1(double);
extern double f2(double);
extern double f3(double);
extern double df1(double);
extern double df2(double);
extern double df3(double);
extern double f4(double);
extern double f5(double);
extern double f6(double);
extern double f7(double);
extern double f8(double);
extern double f9(double);

double root(double func1(double), double func2(double), double a, double b, double eps)
{
    double x = (a + b)/2;
    double g(double y)
    {
        return func1(y) - func2(y);
    }
    if (g(b) < 0)
    {
        double t = b;
        b = a;
        a = t;
    }
    while(fabs(g(x)) > eps)
    {
        iterations++;
        if (g(x)*g(a) < 0)
            b = x;
        else
            a = x;
        x = (a + b)/2;
    }
    return x;
}

double integral(double func(double), double a, double b, double eps)
{
    int n = N_START;
    double h = (b - a)/n;
    double sum2, sum1;
    do
    {
    n *= 2;
    double sum_even = 0;
    double sum_odd = 0;
    for(int i = 1; i < n; i += 1)
        if (i % 2 == 0)
            sum_even += func(a + h*i);
        else
            sum_odd += func(a + h*i);
    sum_even *= 2;
    sum_odd *= 4;
    sum1 = (sum_odd + sum_even + func(a) + func(b))*h/3;
    sum_even = 0;
    sum_odd = 0;
    h = (b - a)/(2*n);
    for(int i = 1; i < 2*n; i += 1)
        if (i % 2 == 0)
            sum_even += func(a + h*i);
        else
            sum_odd += func(a + h*i);
    sum_even *= 2;
    sum_odd *= 4;
    sum2 = (sum_odd + sum_even + func(a) + func(b))*h/3;
    } while(fabs(sum2 - sum1)/CONST_P > eps);
    return sum2;
}

double test_root_eps()
{
    double eps_root = EPS_ROOT_START;
    double s;
    do
    {
    eps_root /= 2;
    double delta_x1 = eps_root/fmin(fabs(df1(X1_MIN) - df2(X1_MIN)),fabs(df1(X1_MAX) - df2(X1_MAX)));
    double delta_x2 = eps_root/fmin(fabs(df1(X2_MIN) - df3(X2_MIN)),fabs(df1(X2_MAX) - df3(X2_MAX)));
    double delta_x3 = eps_root/fmin(fabs(df2(X3_MIN) - df3(X3_MIN)),fabs(df2(X3_MAX) - df3(X3_MAX)));
    double x1 = root(f1, f2, X1_A, X1_B, eps_root);
    double x2 = root(f1, f3, X2_A, X2_B, eps_root);
    double x3 = root(f2, f3, X3_A, X3_B, eps_root);
    double s1, s2, s3, s4, s5, s6;
    s1 = integral(f1, x1, x1 + delta_x1, TEST_EPS_INTEGRAL);
    s2 = integral(f1, x2 - delta_x2, x2, TEST_EPS_INTEGRAL);
    s3 = integral(f2, x1, x1 + delta_x1, TEST_EPS_INTEGRAL);
    s4 = integral(f2, x3 - delta_x3, x3, TEST_EPS_INTEGRAL);
    s5 = integral(f3, x3, x3 + delta_x3, TEST_EPS_INTEGRAL);
    s6 = integral(f3, x2 - delta_x2, x2, TEST_EPS_INTEGRAL);
    s = s1 + s2 + s3 + s4 + s5 + s6;
    } while(s > EPS);
    return eps_root;
}

void test_root()
{
    double x1 = root(f4, f5, TEST_ROOT_X1_A, TEST_ROOT_X1_B, EPS_ROOT);
    double x2 = root(f6, f7, TEST_ROOT_X2_A, TEST_ROOT_X2_B, EPS_ROOT);
    double x3 = root(f8, f9, TEST_ROOT_X3_A, TEST_ROOT_X3_B, EPS_ROOT);
    printf("Root test - sin(x)=cos(x):\n%lf\n", x1);
    printf("Root test - 2^x+2=4^x:\n%lf\n", x2);
    printf("Root test - log2(x-1)=-log2(x+1):\n%lf\n", x3);
}

void test_integral()
{
    double s1 = integral(f4, TEST_INTEGRAL_S1_A, TEST_INTEGRAL_S1_B, EPS_INTEGRAL);
    double s2 = integral(f6, TEST_INTEGRAL_S2_A, TEST_INTEGRAL_S2_B, EPS_INTEGRAL);
    double s3 = integral(df3, TEST_INTEGRAL_S3_A, TEST_INTEGRAL_S3_B, EPS_INTEGRAL);
    printf("Integral test - sin(x), on [0, pi]:\n%lf\n", s1);
    printf("Integral test - 2^x+2, on [0, 2]:\n%lf\n", s2);
    printf("Integral test - 5/x^2, on [2, 4]:\n%lf\n", s3);
}
