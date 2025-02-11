#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <progonka.h>
#define number 20
#define start 1.0
#define finish 2.0
#define counter 10
// double Ai(int i, double *x_points, double step)
// {
//     return ((x_points[i])/(step*step) - (1/(x_points[i])));
// }
// double Bi(int i, double *x_points, double step)
// {
//     return ((x_points[i])/(step*step) + (1/(x_points[i])));
// }
// double Ci(int i, double *x_points, double step)
// {
//     return (-(2*(x_points[i])/(step*step)) - (x_points[i]));
// }
// double Fi(int i, double *x_points, double step)
// {
//     return 0;
// }
double Ai(int i, double *x_points, double step)
{
    return ((x_points[i]) - step);
}
double Bi(int i, double *x_points, double step)
{
    return ((x_points[i]) + step);
}
double Ci(int i, double *x_points, double step)
{
    return (-2 * x_points[i] - x_points[i] * step * step);
}
double Fi(int i, double *x_points, double step)
{
    return 0;
}
double utoch(double x)
{
    return (1.0 / x) * exp(x);
}
double utoch_der(double x)
{
    return (1.0 / x) * exp(x) - (1.0 / (x * x)) * exp(x);
}
void x_filling(double *x_points, int n)
{
    double step = (finish - start) / (n - 1), x = start;
    x_points[0] = start;
    // printf("%.15lf", step);
    for (int i = 1; i < n; ++i)
    {
        x += step;
        x_points[i] = x;
    }
}
void file_write(FILE *file_graph, double *points_x, double *points_y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        fprintf(file_graph, "%.15lf %.15lf\n", points_x[i], points_y[i]);
    }
}
double max_error(double *points_x, double *points_y, int n)
{
    double max = 0;
    for (int i = 0; i < n; ++i)
    {
        if (fabs(points_y[i] - utoch(points_x[i])) > max)
            max = fabs(points_y[i] - utoch(points_x[i]));
    }
    return max;
}
void progonka(double* x_points, double *y_points, int n)
{
    double *k_points = malloc(sizeof(double) * n), *v_points = malloc(sizeof(double) * n);
    double step = (finish - start) / (n - 1);
    x_filling(x_points, n);
    y_points[n - 1] = utoch(x_points[n - 1]);
    k_points[0] = 1.0/(1.0+0.5*step*step);
    v_points[0] = (-step * utoch_der(x_points[0]) + step*step*(utoch_der(x_points[0])))/(1.0 + 0.5*step*step);
      for (int i = 1; i < n; ++i)
    {
        k_points[i] = -(Bi(i, x_points, step) / (Ci(i, x_points, step) + Ai(i, x_points, step) * k_points[i - 1]));
    }
    for (int i = 1; i < n; ++i)
    {
        v_points[i] = -((Ai(i, x_points, step) * v_points[i - 1] - Fi(i, x_points, step)) / (Ci(i, x_points, step) + Ai(i, x_points, step) * k_points[i - 1]));
    }
    for (int i = n - 2; i >= 0; --i)
    {
        y_points[i] = k_points[i] * y_points[i + 1] + v_points[i];
    }
}
int main()
{
    FILE *f = fopen("test.txt", "w");
    // FILE *answer = fopen("answer.txt", "w");
    int n = 2;
    double *x_points = malloc(sizeof(double) * n), *y_points = malloc(sizeof(double) * n), epsilon_old;
    // x_filling(x_points, n);
    // k_points[0] = 1;
    // v_points[0] = -step * utoch_der(x_points[0]);
    // k_points[0] = 0;
    // v_points[0] = utoch(x_points[0]);
    // k_points[0] = 1.0/(1.0+0.5*step*step);
    // v_points[0] = (-step * utoch_der(x_points[0]) + step*step*(utoch_der(x_points[0])))/(1.0 + 0.5*step*step);
    // y_points[n - 1] = utoch(x_points[n - 1]);
    // k_points[0] = (4 + x_points[1] * (-step * step - 2) / (x_points[1] + step)) / (3 - (x_points[1] - step) / (x_points[1] + step));
    // v_points[0] = -(2 * step * utoch_der(x_points[0])) / (3 - (x_points[1] - step) / (x_points[1] + step));
    // progonka(x_points, y_points, n);
    // file_write(f, x_points, y_points, n);
    // printf("\n%.15lf\n", max_error(x_points, y_points, n));
    for (int i = 1; i < counter + 1; ++i)
    {
        if (i != 1)
        {
            x_points = malloc(sizeof(double) * n);
            y_points = malloc(sizeof(double) * n);
            progonka(x_points, y_points, n);
            printf("%.6lf\n", log2(epsilon_old/max_error(x_points, y_points, n)));
            //printf("%.15lf\n", max_error(x_points, y_points, n));
            epsilon_old = max_error(x_points, y_points, n);
            n = n * 2 - 1;
        }
        else
        {
            x_points = malloc(sizeof(double) * n);
            y_points = malloc(sizeof(double) * n);
            progonka(x_points, y_points, n);
            epsilon_old = max_error(x_points, y_points, n);
            n = n * 2 - 1;
        }
    }
    return 0;
}