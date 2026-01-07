#include "butterworth.h"
#include <stdio.h>
#include <time.h>
#include <windows.h>

#define TEST_COUNT 100000 // 测试次数
#define TEST_ANGLE 0x4000  // 测试角度（示例值，Q16格式）

void test_speed(void);

char order = 3;
short bw = 1000;
int b[4] = {0};
int a[4] = {0};

int main()
{
    // printf("%d\n",test());

    // 根据阶数和截止频率计算巴特沃斯滤波器系数b和a
    butter(order, bw, b, a);
    // printf("%d\n",test());
    printf("计算得到的b系数: [");
    for (int i = 0; i <= order; i++)
    {
        printf("%d%s", b[i], (i < order) ? ", " : "]\n");
    }
    printf("计算得到的a系数: [");
    for (int i = 0; i <= order; i++)
    {
        printf("%d%s", a[i], (i < order) ? ", " : "]\n");
    }

    // 性能测试
    test_speed();
    // 导入频率响应数据csv

    // 计算频率响应

    // 将原结果和计算结果输出到新csv
    return 0;
}

void test_speed(void)
{

    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq); // 获取计时器频率
    if (1)
    {
        // 预热（避免首次调用开销）
        butter(order, bw, b, a);
        // 开始计时
        QueryPerformanceCounter(&start);
        // 循环测试
        for (int i = 0; i < TEST_COUNT; i++)
        {
            butter(order, bw, b, a);
        }
        // 结束计时
        QueryPerformanceCounter(&end);
    }
    else
    {
        // 预热（避免首次调用开销）
        fix_tan_pi2_q16(TEST_ANGLE);
        // 开始计时
        QueryPerformanceCounter(&start);
        // 循环测试
        for (int i = 0; i < TEST_COUNT; i++)
        {
            fix_tan_pi2_q16(TEST_ANGLE);
        }
        // 结束计时
        QueryPerformanceCounter(&end);
    }
    // 计算总耗时（秒）
    double total_time = (double)(end.QuadPart - start.QuadPart) / freq.QuadPart;

    // 输出结果（自动选择合适单位）
    if (total_time > 1.0)
    {
        printf("测试 %d 次，总耗时: %.3f 秒\n", TEST_COUNT, total_time);
    }
    else if (total_time > 1e-3)
    {
        printf("测试 %d 次，总耗时: %.3f 毫秒\n", TEST_COUNT, total_time * 1000);
    }
    else
    {
        printf("测试 %d 次，总耗时: %.3f 微秒\n", TEST_COUNT, total_time * 1e6);
    }

    // 计算单次调用耗时（纳秒）
    double avg_time = total_time * 1e9 / TEST_COUNT;
    printf("平均每次调用耗时: %.2f 纳秒\n", avg_time);
}