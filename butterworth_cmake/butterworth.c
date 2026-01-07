/**
 * \file butterworth.c
 * \author changle liu (phantomfancy@outlook.com)
 * \brief Butterworth滤波器实现文件
 * \version v0.5
 * \date 2025-09-29
 * \copyright Copyright (c) 2025, changle liu. All right reserved.
 *
 * \par ChangeLog:
 */
#include "butterworth.h"

#define NULL ((void *)0)
 // q16.16数据类型定义
#define Q16_16_I_MAX 0x7FFF
#define Q16_16_D_MAX 0xFFFF
#define Q16_16_MAX 0x7FFFFFFF
#define Q16_16_ZERO 0x80000000
#define Q16_16_MIN 0xFFFFFFFF

// 数学常数定义（Q16.16 定点格式）
#define PI_Q16_16 0x0003243F               // π ≈ 3.14159
#define PI_HALF_Q16_16 0x0001921F          // π/2 ≈ 1.57079
#define PI_QUARTER_Q16_16 0x0000C910       // π/4 ≈ 0.785398
#define ONE_Q16_16 0x00010000
#define TWO_Q16_16 0x00020000
#define THREE_Q16_16 0x00030000
#define HALF_Q16_16 0x00008000             // 0.5
#define ZERO_Q16_16 0x00000000
#define SQRT2_Q16_16 0x00016ADA            // √2 ≈ 1.41421
#define FIX_WN_THRESHOLD_Q16_16 0x00006666
#define MIN_VAL 0x00000001


#define GET_Q16_16_INT(x) ((short)(x >> 16))
#define GET_Q16_16_DEC(x) ((unsigned short)(x & Q16_16_D_MAX))
#define OPP_Q16_16(x) (-(x))
#define ABS_Q16_16(x) (((x) < 0) ? (-(x)) : (x))
static inline q16_16_t MUL_Q16_16(q16_16_t a, q16_16_t b)
{
    if (a == 0x00010000) return b;
    if (b == 0x00010000) return a;
    return (q16_16_t)(((long long)a * (long long)b) >> 16);
}
static inline q16_16_t DIV_Q16_16(q16_16_t a, q16_16_t b)
{
    if (a == b) return 0x00010000;
    if (b == 0x00010000) return a;
    return (q16_16_t)(((long long)a << 16) / (b == 0 ? 1 : b));
}
#define TRANS_INT_TO_Q16_16(x) ((int)((x >= 0) ? (x << 16) : ((-(x)) << 16)))
#define TRANS_Q16_16_TO_INT(x) ((int)(GET_Q16_16_INT(x) * 10000 + ((GET_Q16_16_DEC(x) * 10000 + 0x8000) >> 16)))

// int invoke_count = 0;

q16_16_t fix_tan_pi2_q16(q16_16_t Wn_q16)
{
    // invoke_count++;
    if (Wn_q16 < 0x00002000) {  // Wn < 0.1,tanx ≈ x
        return MUL_Q16_16(PI_HALF_Q16_16, Wn_q16);
    }
    else if (Wn_q16 <= FIX_WN_THRESHOLD_Q16_16) // 在泰勒级数的误差小的区间，使用霍纳法则分层计算
    {
        q16_16_t theta = MUL_Q16_16(PI_HALF_Q16_16, Wn_q16); // θ = (π/2) * Wn
        q16_16_t t2 = MUL_Q16_16(theta, theta);
        return MUL_Q16_16(theta, (ONE_Q16_16 + MUL_Q16_16(t2, (0x00005555 + MUL_Q16_16(t2, (0x00002222 + MUL_Q16_16(t2, (0x00000DD1 + MUL_Q16_16(t2, 0x00000599)))))))));
    }
    else // 否则递归计算tan(x/2) = 2tan(x/2)/(1-tan^2(x/2))
    {
        q16_16_t saturate_half_tan_pi2 = fix_tan_pi2_q16(MUL_Q16_16(HALF_Q16_16, Wn_q16));
        return DIV_Q16_16(MUL_Q16_16(TWO_Q16_16, saturate_half_tan_pi2), (ONE_Q16_16 - MUL_Q16_16(saturate_half_tan_pi2, saturate_half_tan_pi2)));
    }
}

char butter(unsigned char order, short bw, int* b, int* a)
{
    // 输入参数检查
    if ((order < 2) || (order > BW_MAX_ORDER) || (b == NULL) || (a == NULL) || (bw <= 0) || (bw > LOOP_FREQ_HZ / 2))
    {
        return -1;
    }

    // 步骤1: 计算预畸变截止频率 wc = 2 * tan(pi/2 * Wn),以及其他参数
    q16_16_t wc = MUL_Q16_16(TWO_Q16_16, fix_tan_pi2_q16((q16_16_t)((((long long)bw * 2 * 0x10000LL) / LOOP_FREQ_HZ))));
    q16_16_t proto[BW_MAX_ORDER + 1] = { 0 }; // butterworth prototype
    q16_16_t K = TWO_Q16_16;
    q16_16_t K2 = MUL_Q16_16(K, K);
    if (order == 2)
    {
        // 步骤2: 模拟域传递函数
        // 2阶巴特沃斯原型: H(s) = 1 / (s² + sqrt(2)s + 1)
        proto[0] = MUL_Q16_16(wc, wc);
        proto[1] = MUL_Q16_16(SQRT2_Q16_16, wc);
        proto[2] = ONE_Q16_16;

        // 步骤3: 双线性变换
        q16_16_t inv_K2 = DIV_Q16_16(ONE_Q16_16, K2);
        a[0] = proto[0] + MUL_Q16_16(proto[1], K) + K2;
        a[1] = MUL_Q16_16(TWO_Q16_16, proto[0]) - MUL_Q16_16(TWO_Q16_16, K2);
        a[2] = proto[0] - MUL_Q16_16(proto[1], K) + K2;

        b[0] = inv_K2;
        b[1] = MUL_Q16_16(TWO_Q16_16, inv_K2);
        b[2] = inv_K2;
    }
    else if (order == 3)
    {
        // 步骤2: 模拟域传递函数
        // 3阶巴特沃斯原型: H(s) = 1 / (s³ + 2s² + 2s + 1)
        proto[0] = MUL_Q16_16(wc, wc);
        proto[1] = MUL_Q16_16(TWO_Q16_16, proto[0]);
        proto[0] = MUL_Q16_16(proto[0], wc);
        proto[2] = MUL_Q16_16(TWO_Q16_16, wc);
        proto[3] = ONE_Q16_16;

        // 步骤3: 双线性变换
        q16_16_t K3 = MUL_Q16_16(K2, K);
        q16_16_t inv_K3 = DIV_Q16_16(ONE_Q16_16, K3);
        a[0] = proto[0] + MUL_Q16_16(proto[1], K) + MUL_Q16_16(proto[2], K2) + K3;
        a[1] = MUL_Q16_16(THREE_Q16_16, proto[0]) + MUL_Q16_16(proto[1], K) - MUL_Q16_16(proto[2], K2) - MUL_Q16_16(THREE_Q16_16, K3);
        a[2] = MUL_Q16_16(THREE_Q16_16, proto[0]) - MUL_Q16_16(proto[1], K) - MUL_Q16_16(proto[2], K2) + MUL_Q16_16(THREE_Q16_16, K3);
        a[3] = proto[0] - MUL_Q16_16(proto[1], K) + MUL_Q16_16(proto[2], K2) - K3;

        b[0] = inv_K3;
        b[1] = MUL_Q16_16(THREE_Q16_16, inv_K3);
        b[2] = MUL_Q16_16(THREE_Q16_16, inv_K3);
        b[3] = inv_K3;
    }

    // 步骤4: DC增益归一化
    q16_16_t dc_num = ZERO_Q16_16, dc_den = ZERO_Q16_16;
    for (int i = 0; i <= order; ++i)
    {
        dc_num += b[i];
        dc_den += a[i];
    }

    // 步骤5: 系数归一化（确保a[0] = 1）和符号修正
    q16_16_t norm = ONE_Q16_16;
    if ((dc_num > MIN_VAL) && (dc_den > MIN_VAL))
    {
        norm = DIV_Q16_16(dc_den, dc_num);
    }
    q16_16_t sign_corr = (a[0] < ZERO_Q16_16) ? OPP_Q16_16(ONE_Q16_16) : ONE_Q16_16;
    for (int i = order; i >= 0; --i)
    {
        b[i] = MUL_Q16_16(DIV_Q16_16(MUL_Q16_16(b[i], norm), a[0]), sign_corr);
        a[i] = MUL_Q16_16(DIV_Q16_16(a[i], a[0]), sign_corr);
        b[i] = TRANS_Q16_16_TO_INT(b[i]);
        a[i] = TRANS_Q16_16_TO_INT(a[i]);
    }
    return 0;
}