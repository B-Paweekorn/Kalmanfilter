/*
 * kalman.c
 *
 *  Created on: 17 May 2023
 *      Author: Paweekorn
 */

#include "kalman.h"
#include "arm_math.h"
float Kalman_Speed = 0;

// Define initial state of X_pk and P_pk
float32_t X_k[4] = {0.0f, 0.0f, 0.0f, 0.0f};
arm_matrix_instance_f32 X_k_matrix;

float32_t P_k[16] = {0.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f};
arm_matrix_instance_f32 P_k_matrix;

// System matrices
//float32_t A[16] = {1.0f, 0.0091f, -0.0516f, 0.0047f,
//                   0.0f, 0.7765f, -9.8835f, 0.6533f,
//                   0.0f, 0.0f   ,   1.0f  , 0.0f,
//                   0.0f,-0.1651f, 1.2950f, 0.0045f};
//arm_matrix_instance_f32 A_matrix;
//
//float32_t A_transpose[16] = {1.0f,    0.0f,    0.0f,   0.0f,
//	  	  	  	  	  	 0.0091f, 0.7765f, 0.0f, -0.1651f,
//	  	  	  	  	  	-0.0516f,-9.8835f, 1.0f,  1.2950f,
//	  	  	  	  	  	 0.0047f, 0.6533f,0.0f,  0.0045f};
//arm_matrix_instance_f32 A_transpose_matrix;

//float32_t A[16] = {1.0f, 1.9998e-4f, -2.1683e-5f, 3.6545e-6f,
//                   0.0f, 0.9998f, -0.2168f, 0.0363f,
//                   0.0f, 0.0f   ,   1.0f  , 0.0f,
//                   0.0f,-0.0092f, 0.0010f, 0.9569f};
//arm_matrix_instance_f32 A_matrix;
//
//float32_t A_transpose[16] = {1.0f,    0.0f,    0.0f,   0.0f,
//	  	  	  	  	  	 1.9998e-4f, 0.9998f, 0.0f, -0.0092f,
//	  	  	  	  	  	-2.1683e-5f,-0.2168f, 1.0f,  0.0010f,
//	  	  	  	  	  	 3.6545e-6f, 0.0363f, 0.0f,  0.9569f};
//arm_matrix_instance_f32 A_transpose_matrix;

float32_t A[16] = {1.0f, 0.000996946187806927f,-0.000781560031349208f, 8.37430695770173e-05f,
                   0.0f, 0.992523116469047f,   -1.56115964889452f, 0.160829336670584f,
                   0.0f, 0.0f   ,   1.0f  , 0.0f,
                   0.0f,-0.0454381199874783f,  0.0370492618434313f, 0.782610780493230f};
arm_matrix_instance_f32 A_matrix;

float32_t A_transpose[16] = {1.0f,    	 0.0f,    0.0f,   0.0f,
		0.000996946187806927f, 0.992523116469047f, 0.0f, -0.0454381199874783f,
		-0.000781560031349208f,-1.56115964889452f, 1.0f,  0.0370492618434313f,
		8.37430695770173e-05f, 0.160829336670584f, 0.0f,  0.782610780493230f};
arm_matrix_instance_f32 A_transpose_matrix;


float32_t eye[16] = {1.0f, 0.0f, 0.0f, 0.0f,
		  	  	 0.0f, 1.0f, 0.0f, 0.0f,
				 0.0f, 0.0f, 1.0f, 0.0f,
				 0.0f, 0.0f, 0.0f, 1.0f,};
arm_matrix_instance_f32 eye_matrix;

//float32_t B[4] = {4.7925e-8f,
//				  7.162e-4f,
//				  0.0f,
//				  0.0384f};
//arm_matrix_instance_f32 B_matrix;

float32_t B[4] = {6.11011237621270e-06f,
			      0.0179660664417631f,
				  0.0f,
				  0.190433717271840f};
arm_matrix_instance_f32 B_matrix;

float32_t C[4] = {1.0f, 0.0f, 0.0f, 0.0f};
arm_matrix_instance_f32 C_matrix;

float32_t C_transpose[4] = {1.0f,
	  	  	  	  	  	    0.0f,
	  	  	  	  	  	    0.0f,
  						    0.0f};
arm_matrix_instance_f32 C_transpose_matrix;

float32_t G[4] = {0.0f,
				  1.0f,
				  0.0f,
				  0.0f};
arm_matrix_instance_f32 G_matrix;

float32_t G_transpose[4] = {0.0f,1.0f,0.0f,0.0f};
arm_matrix_instance_f32 G_transpose_matrix;

float32_t Es_velocity[1] = {0.0f};
arm_matrix_instance_f32 Output_matrix;

arm_matrix_instance_f32 GGT_matrix;
float32_t GGT[16];

arm_matrix_instance_f32 GQGT_matrix;
float32_t GQGT[16];
//------------------------------------------- for equa ------------------------------------
// Compute Xk = Ax + Bu
arm_matrix_instance_f32 Bu_matrix;
float32_t Bu_data[4];
arm_matrix_instance_f32 Ax_matrix;
float32_t Ax_data[4];

// Compute (C * P_k * C^T + R)
arm_matrix_instance_f32 CP_matrix;
float32_t CP[4];
arm_matrix_instance_f32 CPCT_matrix;
float32_t CPCT[1];
arm_matrix_instance_f32 CPCTR_matrix;
float32_t CPCTR[1];

// Compute Kalman Gain: K = P_k * C^T * inv(C * P_k * C^T + R)
arm_matrix_instance_f32 K_matrix;
float32_t K[4];
arm_matrix_instance_f32 PCT_matrix;
float32_t PCT[4];

// Compute inverse of (C * P_k * C^T + R)
arm_matrix_instance_f32 CPCTRinv_matrix;
float32_t CPCTRinv[1];

// Computation of the estimated state
arm_matrix_instance_f32 Cx_matrix;
float32_t Cx[1];
arm_matrix_instance_f32 yCx_matrix;
float32_t yCx[1];
arm_matrix_instance_f32 KyCx_matrix;
float32_t KyCx[4];


float32_t Q = 1.0;
arm_matrix_instance_f32 R_matrix;
float32_t R[1] = {1.0f};

arm_matrix_instance_f32 Z_matrix;
float32_t Z[1];

arm_matrix_instance_f32 Velocity_matrix;

volatile arm_status Calst;

float checkVal;

float SteadyStateKalmanFilter(float32_t Vin,float32_t Velocity){
	  arm_mat_init_f32(&Velocity_matrix, 1, 1,(float32_t*) &Velocity);
	  // Compute Xk = Ax + Bu
	  arm_mat_scale_f32(&B_matrix, Vin, &Bu_matrix); 		   				// Bu
	  arm_mat_mult_f32(&A_matrix, &X_k_matrix, &Ax_matrix);  		   		// Ax
	  arm_mat_add_f32(&Ax_matrix, &Bu_matrix, &X_k_matrix); 		   		// Xk = Ax + Bu

	  // Compute (A * P_pk * A^T + G * Q * G^T)
	  arm_mat_mult_f32(&A_matrix, &P_k_matrix, &P_k_matrix);  		   		// Pk = A * P_pk
	  arm_mat_mult_f32(&P_k_matrix, &A_transpose_matrix, &P_k_matrix); 		// Pk = A * P_pk * A^T
	  arm_mat_mult_f32(&G_matrix, &G_transpose_matrix, &GGT_matrix);        // G * G^T
	  arm_mat_scale_f32(&GGT_matrix, Q, &GQGT_matrix); 				   	   	// G * Q
	  arm_mat_add_f32(&P_k_matrix, &GQGT_matrix, &P_k_matrix); 	       		// A * P_pk * A^T + G * Q * G^T

	  // Compute (C * P_k * C^T + R)
	  arm_mat_mult_f32(&C_matrix, &P_k_matrix, &CP_matrix);			     // C * Pk
	  arm_mat_mult_f32(&CP_matrix, &C_transpose_matrix, &CPCT_matrix);   // C * Pk * C^T
	  arm_mat_add_f32(&CPCT_matrix, &R_matrix, &CPCTR_matrix);			 // C * P_k * C^T + R

	  // Compute inverse of (C * P_k * C^T + R)
	  arm_mat_inverse_f32(&CPCTR_matrix, &CPCTRinv_matrix);					 // inverse of (C * P_k * C^T + R)

	  // Compute Kalman Gain: K = P_k * C^T * inv(C * P_k * C^T + R)
	  arm_mat_mult_f32(&P_k_matrix, &C_transpose_matrix, &PCT_matrix); 		 // P_k * C^T
	  arm_mat_mult_f32(&PCT_matrix, &CPCTRinv_matrix, &K_matrix);  			 // P_k * C^T * inv(C * P_k * C^T + R)

	  // Computation of the estimated state
	  arm_mat_mult_f32(&C_matrix, &X_k_matrix, &Cx_matrix);				 // C * X_k
	  arm_mat_sub_f32(&Velocity_matrix,  &Cx_matrix, &yCx_matrix);			  // y - ( C * X_k )
	  arm_mat_mult_f32(&K_matrix, &yCx_matrix, &KyCx_matrix);		     // K( y - ( C * X_k ) )
	  arm_mat_add_f32(&X_k_matrix, &KyCx_matrix, &X_k_matrix);		 	 // X_k + K( y - ( C * X_k ) )

	  // Computation of the estimated output
	  arm_mat_mult_f32(&C_matrix, &X_k_matrix, &Output_matrix);

	  // Computation of the state covariance error
	  arm_matrix_instance_f32 temp_matrix4;
	  float32_t temp_data4[16];
	  arm_mat_init_f32(&temp_matrix4, 4, 4,(float32_t*) &temp_data4);

	  arm_mat_mult_f32(&K_matrix, &C_matrix, &temp_matrix4);				// K * C
	  arm_mat_sub_f32(&eye_matrix, &temp_matrix4, &temp_matrix4);			// (I - (K * C))
	  arm_mat_mult_f32(&temp_matrix4, &P_k_matrix, &P_k_matrix);			// (I - (K * C)) * P_k
	  Kalman_Speed = X_k[1];
	  return  Kalman_Speed;
}

void Kalman_Start(){
	arm_mat_init_f32(&X_k_matrix, 4, 1,(float32_t*) &X_k);
	arm_mat_init_f32(&P_k_matrix, 4, 4,(float32_t*) &P_k);

	arm_mat_init_f32(&A_matrix, 4, 4,(float32_t*) &A);
	arm_mat_init_f32(&B_matrix, 4, 1,(float32_t*) &B);
	arm_mat_init_f32(&C_matrix, 1, 4,(float32_t*) &C);
	arm_mat_init_f32(&G_matrix, 4, 1,(float32_t*) &G);

	arm_mat_init_f32(&A_transpose_matrix, 4, 4,(float32_t*) &A_transpose);
	arm_mat_init_f32(&C_transpose_matrix, 4, 1,(float32_t*) &C_transpose);
	arm_mat_init_f32(&G_transpose_matrix, 1, 4,(float32_t*) &G_transpose);

	arm_mat_init_f32(&GGT_matrix, 4, 4,(float32_t*) &GGT);
	arm_mat_init_f32(&GQGT_matrix, 4, 4,(float32_t*) &GQGT);

	// Compute Xk = Ax + Bu
	arm_mat_init_f32(&Bu_matrix, 4, 1,(float32_t*) &Bu_data);
	arm_mat_init_f32(&Ax_matrix, 4, 1,(float32_t*) &Ax_data);

	// Compute (C * P_k * C^T + R)
	arm_mat_init_f32(&CP_matrix, 1, 4,(float32_t*) &CP);
	arm_mat_init_f32(&CPCT_matrix, 1, 1,(float32_t*) &CPCT);
	arm_mat_init_f32(&CPCTR_matrix, 1, 1,(float32_t*) &CPCTR);

	// Compute Kalman Gain: K = P_k * C^T * inv(C * P_k * C^T + R)
	arm_mat_init_f32(&K_matrix, 4, 1,(float32_t*) &K);
	arm_mat_init_f32(&PCT_matrix, 4, 1,(float32_t*) &PCT);

	// Compute inverse of (C * P_k * C^T + R)
	arm_mat_init_f32(&CPCTRinv_matrix, 1, 1,(float32_t*) &CPCTRinv);

	// Computation of the estimated state
	arm_mat_init_f32(&Cx_matrix, 1, 1,(float32_t*) &Cx);
	arm_mat_init_f32(&yCx_matrix, 1, 1,(float32_t*) &yCx);
	arm_mat_init_f32(&KyCx_matrix, 4, 1,(float32_t*) &KyCx);

	arm_mat_init_f32(&Output_matrix, 1, 1,(float32_t*) &Es_velocity);

	arm_mat_init_f32(&eye_matrix, 4, 4,(float32_t*) &eye);

	arm_mat_init_f32(&R_matrix, 1, 1,(float32_t*) &R);
	arm_mat_init_f32(&Z_matrix, 1, 1,(float32_t*) &Z);
}



