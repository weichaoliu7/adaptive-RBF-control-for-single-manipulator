#include <stdio.h>
#include <stdlib.h>
#include "sine.h"
#include "cosine.h"

#define PI 3.14159

// reference: [1]Liu JinKun. Robot Control System Design and MATLAB Simulation[M]. Tsinghua University Press, 2008.
// [2]Feng G. A compensating scheme for robot tracking based on neural networks[J]. Robotics and Autonomous Systems, 1995, 15(3): 199-206.

void matrix_Multi(double C[][2], double A[][2], double B[][2], int rows1, int cols1, int cols2){
    for (int i = 0; i < rows1; i++){
        for (int j = 0; j < cols2; j++){
            C[i][j] = 0;
            for (int k = 0; k < cols1; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main()
{

    double Ts = 0.005; // sampling period
    double t0 = 0.0;
    double t1 = 30.0;
    double qd_amp = 1.0;
    double dqd_amp = 1.0;
    double ddqd_amp = -1.0;
    double d_amp = 1.3;

    Data qd, ddqd, d;
    Data1 dqd;
    sine(&qd, 1, Ts, t0, t1, qd_amp);     // desired angular displacement
    cosine(&dqd, 1, Ts, t0, t1, dqd_amp); // desired angular velocity
    sine(&ddqd, 1, Ts, t0, t1, ddqd_amp); // desired angular acceleration

    int ARRAY_SIZE = (t1 - t0) / Ts;

    int H = 5, IN = 2; // hidden layer cell and input layer cell's number
    int M1 = 2;        // adaptive law selection
    int M = 3;         // control scheme selection for precise model
    double gamma = 150;
    double kp = 20, kv = 10; // proportional gain and differential gain
    double ctrl_u1[ARRAY_SIZE], ctrl_u2[ARRAY_SIZE], ctrl_u3[ARRAY_SIZE], ctrl_u4[ARRAY_SIZE], plant_u[ARRAY_SIZE];
    double q[ARRAY_SIZE], dq[ARRAY_SIZE], ddq[ARRAY_SIZE];
    double f[ARRAY_SIZE], fn[ARRAY_SIZE], tol1[ARRAY_SIZE], tol2[ARRAY_SIZE], tol[ARRAY_SIZE];
    double e, de, x[2], h[H], c[IN][H], b[H], theta[H], d_theta[H];
    double h_x[H][IN], h_x_P[H][IN], h_x_P_B[H];

    for (int i = 0; i < ARRAY_SIZE; i++){

        if (i == 0){
            for (int j = 0; j < H; j++){
                theta[j] = 0.1;
            }
            q[i] = 0.15;
            dq[i] = 0.0;
            ddq[i] = 0.0;
        }

        // controller
        ctrl_u1[i] = qd.y[i];        // system input
        ctrl_u2[i] = q[i];           // angular displacement
        ctrl_u3[i] = dq[i];          // angular velocity
        ctrl_u4[i] = ddq[i];         // angular acceleration
        e = ctrl_u2[i] - ctrl_u1[i]; // angular displacement error
        de = ctrl_u3[i] - dqd.y[i];  // angular velocity error
        x[0] = e;                    // state variable
        x[1] = de;

        for (int j = 0; j < H; j++){
            for (int k = 0; k < IN; k++){
                c[j][k] = 0.60;
            }
        }
        for (int j = 0; j < H; j++){
            b[j] = 3.0;
        }

        for (int j = 0; j < H; j++){
            double sum = 0.0;
            for (int k = 0; k < IN; k++){
                sum += pow(x[k] - c[j][k], 2);
            }
            h[j] = exp(-sum / (2 * b[j] * b[j])); // RBF function's output
        }

        double A[2][2] = {{0, 1}, {-kp, -kv}};
        double B[2] = {0, 1};
        double Q[2][2] = {{50, 0}, {0, 50}};
        double P[2][2] = {{65, 1.25}, {1.25, 2.625}}; // use matlab get the solution of Lyapunov equation: P=lyap(A',Q);

        for (int j = 0; j < H; j++){
            for (int k = 0; k < IN; k++){
                h_x[j][k] = h[j] * x[k];
            }
        }

        matrix_Multi(h_x_P, h_x, P, H, IN, 2);
        for (int j = 0; j < H; j++){
            h_x_P_B[j] = 0.0;
            for (int k = 0; k < IN; k++){
                h_x_P_B[j] += h_x_P[j][k] * B[k];
            }
        }

        if (M1 == 1){ // Adaptive Law
            for (int j = 0; j < H; j++){
                d_theta[j] = gamma * h_x_P_B[j]; // derivative of weight estimate
            }
        }
        else if (M1 == 2){ // Adaptive Law with UUB
            double k1 = 0.001;
            for (int j = 0; j < H; j++){
                d_theta[j] = gamma * h_x_P_B[j] + k1 * gamma * sqrt(pow(x[0], 2) + pow(x[1], 2)) * theta[j];
            }
        }

        for (int j = 0; j < H; j++){
            theta[j] = theta[j] + d_theta[j] * Ts;
        }

        double g = 9.8;                          // gravitational acceleration
        double m = 1, l = 0.25;                  // mass and length of link
        double D0 = 4.0 / 3.0 * m * pow(l, 2);   // inertia matrix of nominal model
        double C0 = 2.0;                         // coriolis term of the nominal model
        double G0 = m * g * l * cos(ctrl_u2[i]); // gravity term of the nominal model

        double d_D0 = 0.2 * D0;
        double d_C0 = 0.2 * C0;
        double d_G0 = 0.2 * G0;

        double D = D0 - d_D0; // inertia matrix of precise model
        double C = C0 - d_C0; // coriolis term of precise model
        double G = G0 - d_G0; // gravity term of precise model

        sine(&d, 0.5 * PI, Ts, t0, t1, d_amp); // external disturbance

        if (M == 1){                                                                         // control for precise model
            tol1[i] = D0 * (ddqd.y[i] - kv * de - kp * e) + C0 * ctrl_u3[i] + G0; // computed torque controller
            tol2[i] = 0;
            tol[i] = tol1[i] + tol2[i];
        }
        else if (M == 2){                                                                            // control with precise nonlinear compensation
            f[i] = 1 / D0 * (d_D0 * ctrl_u4[i] + d_C0 * ctrl_u3[i] + d_G0 + d.y[i]); // model uncertainty term
            tol1[i] = D0 * (ddqd.y[i] - kv * de - kp * e) + C0 * ctrl_u3[i] + G0;
            tol2[i] = -D0 * f[i];
            tol[i] = tol1[i] + tol2[i]; // improved computed torque controller
        }
        else if (M == 3){ // control with neural compensation
            tol1[i] = D0 * (ddqd.y[i] - kv * de - kp * e) + C0 * ctrl_u3[i] + G0;
            f[i] = 1 / D0 * (d_D0 * ctrl_u4[i] + d_C0 * ctrl_u3[i] + d_G0 + d.y[i]);
            fn[i] = 0.0;
            for (int j = 0; j < H; j++){
                fn[i] += theta[j] * h[j]; // NN's output = transposition of weight estimate * RBF function's output
            }
            tol2[i] = -D0 * fn[i];
            tol[i] = tol1[i] + tol2[i];
        }

        // plant
        int P1 = 2;
        if (P1 == 1){
            ddq[i] = 1 / D0 * (-C0 * dq[i] - G0 + tol[i]);
        }
        else if (P1 == 2){
            ddq[i] = 1 / D * (-C * dq[i] - G + tol[i] + d.y[i]);
        }

        dq[i + 1] = dq[i] + ddq[i] * Ts;
        q[i + 1] = q[i] + dq[i] * Ts;
        ddq[i + 1] = ddq[i];
    }

    FILE *fp1 = fopen("qd.txt", "w");
    FILE *fp2 = fopen("q.txt", "w");
    FILE *fp3 = fopen("tol1.txt", "w");
    FILE *fp4 = fopen("tol2.txt", "w");
    FILE *fp5 = fopen("tol.txt", "w");
    FILE *fp6 = fopen("f.txt", "w");
    FILE *fp7 = fopen("fn.txt", "w");

    if (fp1 != NULL && fp2 != NULL && fp3 != NULL && fp4 != NULL && fp5 != NULL && fp6 != NULL && fp7 != NULL){
        for (int i = 0; i < ARRAY_SIZE; i++){
            fprintf(fp1, "%lf\n", qd.y[i]);
            fprintf(fp2, "%lf\n", q[i]);
            fprintf(fp3, "%lf\n", tol1[i]);
            fprintf(fp4, "%lf\n", tol2[i]);
            fprintf(fp5, "%lf\n", tol[i]);
            fprintf(fp6, "%lf\n", f[i]);
            fprintf(fp7, "%lf\n", fn[i]);
        }
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
        fclose(fp5);
        fclose(fp6);
        fclose(fp7);
    }
    else{
        printf("unable to open file\n");
    }

    return 0;
}
