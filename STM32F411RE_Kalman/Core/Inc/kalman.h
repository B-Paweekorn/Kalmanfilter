/*
 * kalman.h
 *
 *  Created on: May 16, 2023
 *      Author: pawee
 */

#ifndef INC_KALMAN_H_
#define INC_KALMAN_H_

#include "main.h"
#include "arm_math.h"

float SteadyStateKalmanFilter(float32_t Vin, float32_t Velocity);
void Kalman_Start();

#endif /* INC_KALMAN_H_ */
