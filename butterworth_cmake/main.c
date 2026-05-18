#include "butterworth.h"
#include <stdio.h>
#include <time.h>
#include <windows.h>

#define TEST_COUNT 100000 // 测试次数

void test_butter_result(void);
void test_speed(void);
static void print_coefficients(const char* method_name, int* b_coeff, int* a_coeff);
static void print_diff(const char* coeff_name, int* lhs, int* rhs);
static double benchmark_butter_method(butter_tan_method_t tan_method);
static void print_time_result(const char* method_name, double total_time);

char order = 3;
short bw = 1000;
int b[4] = {0};
int a[4] = {0};

int main()
{
    // printf("%d\n",test());

    // 根据阶数和截止频率计算巴特沃斯滤波器系数b和a
    test_butter_result();

    // 性能测试
    test_speed();

    return 0;
}

void test_butter_result(void)
{
    int b_taylor[4] = {0};
    int a_taylor[4] = {0};
    int b_cordic[4] = {0};
    int a_cordic[4] = {0};

    char ret_taylor = butter_with_tan_method(order, bw, b_taylor, a_taylor, BUTTER_TAN_METHOD_TAYLOR);
    char ret_cordic = butter_with_tan_method(order, bw, b_cordic, a_cordic, BUTTER_TAN_METHOD_CORDIC);

    if ((ret_taylor != 0) || (ret_cordic != 0))
    {
        printf("系数计算失败: Taylor=%d, CORDIC=%d\n", ret_taylor, ret_cordic);
        return;
    }

    print_coefficients("Taylor", b_taylor, a_taylor);
    print_coefficients("CORDIC", b_cordic, a_cordic);
    print_diff("b", b_cordic, b_taylor);
    print_diff("a", a_cordic, a_taylor);
}

static void print_coefficients(const char* method_name, int* b_coeff, int* a_coeff)
{
    printf("%s方法计算得到的b系数: [", method_name);
    for (int i = 0; i <= order; i++)
    {
        printf("%d%s", b_coeff[i], (i < order) ? ", " : "]\n");
    }
    printf("%s方法计算得到的a系数: [", method_name);
    for (int i = 0; i <= order; i++)
    {
        printf("%d%s", a_coeff[i], (i < order) ? ", " : "]\n");
    }
}

static void print_diff(const char* coeff_name, int* lhs, int* rhs)
{
    printf("CORDIC - Taylor 的%s系数差值: [", coeff_name);
    for (int i = 0; i <= order; i++)
    {
        printf("%d%s", lhs[i] - rhs[i], (i < order) ? ", " : "]\n");
    }
}

void test_speed(void)
{
    double taylor_time = benchmark_butter_method(BUTTER_TAN_METHOD_TAYLOR);
    double cordic_time = benchmark_butter_method(BUTTER_TAN_METHOD_CORDIC);

    print_time_result("Taylor", taylor_time);
    print_time_result("CORDIC", cordic_time);
    if (taylor_time > 0.0)
    {
        printf("CORDIC / Taylor 耗时倍率: %.2f\n", cordic_time / taylor_time);
    }
}

static double benchmark_butter_method(butter_tan_method_t tan_method)
{
    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq); // 获取计时器频率

    // 预热（避免首次调用开销）
    butter_with_tan_method(order, bw, b, a, tan_method);
    // 开始计时
    QueryPerformanceCounter(&start);
    // 循环测试
    for (int i = 0; i < TEST_COUNT; i++)
    {
        butter_with_tan_method(order, bw, b, a, tan_method);
    }
    // 结束计时
    QueryPerformanceCounter(&end);

    // 计算总耗时（秒）
    return (double)(end.QuadPart - start.QuadPart) / freq.QuadPart;
}

static void print_time_result(const char* method_name, double total_time)
{
    printf("%s方法测试 %d 次，总耗时: ", method_name, TEST_COUNT);
    // 输出结果（自动选择合适单位）
    if (total_time > 1.0)
    {
        printf("%.3f 秒\n", total_time);
    }
    else if (total_time > 1e-3)
    {
        printf("%.3f 毫秒\n", total_time * 1000);
    }
    else
    {
        printf("%.3f 微秒\n", total_time * 1e6);
    }

    // 计算单次调用耗时（纳秒）
    double avg_time = total_time * 1e9 / TEST_COUNT;
    printf("%s方法平均每次调用耗时: %.2f 纳秒\n", method_name, avg_time);
}