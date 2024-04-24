# Example usage

in main.c

1. /* USER CODE BEGIN PV */

```c
KalmanFilter filterA;
```

2. /* USER CODE BEGIN 2 */
```c
Kalman_Start(&filterA);
```

3. In the loop // struct, voltage(V), motor_speed(rad/s)

```c
SteadyStateKalmanFilter(&filterA, volt[i], motor_speed[i] * 2 * 3.14);
```
