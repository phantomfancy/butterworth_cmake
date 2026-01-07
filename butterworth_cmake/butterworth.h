/**
 * \file butterworth.h
 * \author changle liu (phantomfancy@outlook.com)
 * \brief Butterworth滤波器头文件
 * \version v0.5
 * \date 2025-09-29
 * \copyright Copyright (c) 2025, changle liu. All right reserved.
 *
 * \par ChangeLog:
 */
#pragma once

 /**
  * \brief q格式定义.
  * q[num1]_[num2]_t 代表num1.num2的q格式数,其中num1为整数部分位数,num2为小数部分位数.
  * 整数部分为有符号数,小数部分为无符号数.
  */
typedef int q16_16_t;

#define BW_MAX_ORDER 3     // 最大滤波器阶数
#define LOOP_FREQ_HZ 10000 // 采样频率(Hz)

/**
 * \brief 计算[-π/2, π/2]范围内的正切值（Q16.16定点格式）
 * \param Wn_q16 输入角度（Q16.16格式）
 * \return 正切值（Q16.16格式）
 */
q16_16_t fix_tan_pi2_q16(q16_16_t Wn_q16);

/**
 * \brief  设计Butterworth数字滤波器
 * @param  order: 滤波器阶数 (1~BW_MAX_ORDER)
 * @param  bw: 截止频率(0 ~ LOOP_FREQ_HZ/2 Hz)
 * @param  b: 输出分子系数数组
 * @param  a: 输出分母系数数组
 * \retval 成功返回0，失败返回错误码
 */
char butter(unsigned char order, short bw, int* b, int* a);