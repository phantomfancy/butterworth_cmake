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
#pragma once

 /**
  * \brief q格式定义.
  * qnum1_num2_t 代表num1.num2的q格式数,其中num1为整数部分位数,num2为小数部分位数.
  * 整数部分为有符号数,小数部分为无符号数.使用uqnum1_num2_t时,整数部分为无符号数.
  */
  // 0.16格式无符号定点数
  // typedef unsigned short q0_16_t;
  // // 0.32格式无符号定点数
  // typedef unsigned int q0_32_t;
  // // 1.15格式有符号定点数
  // typedef short q15_t;
  // // 1.31格式有符号定点数
  // typedef int q31_t;
  // // 8.8格式有符号定点数
  // typedef short q8_8_t;
  // // 16.16格式有符号定点数
typedef int q16_16_t;
// // 32.32格式有符号定点数
// typedef long long q32_32_t;

#define BW_MAX_ORDER 3     // 最大滤波器阶数
#define LOOP_FREQ_HZ 10000 // 采样频率(Hz)

/**
 * \brief 计算[-π/2, π/2]范围内的正切值（Q16.16定点格式）
 * \param Wn_q16 输入角度（Q16.16格式）
 * \return 正切值（Q16.16格式）
 */
q16_16_t fix_tan_pi2_q16(q16_16_t Wn_q16);

/**
 * \brief  使用芯片上下文参数设计Butterworth滤波器系数
 * @param  order: 滤波器阶数 (1~BW_MAX_ORDER)
 * @param  bw: 截止频率(0 ~ LOOP_FREQ_HZ/2 Hz)
 * @param  b: 输出分子系数数组
 * @param  a: 输出分母系数数组
 * \retval 成功返回0，失败返回错误码
 */
char butter(unsigned char order, short bw, int* b, int* a);

/**
 * \brief  for pc test only!
 * \retval int 
 */
int test(void);