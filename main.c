#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <progonka.h>
// Добавить вычисление первого шага
// Проверить прогонку
#define counter 10
double Ai(double time_step, double space_step)
{
    return -time_step / space_step;
}
double Bi(double time_step, double space_step)
{
    return time_step / space_step;
}
double Ci(double time_step, double space_step)
{
    return 3.0;
}
double Fi(int n, int j, double **u)
{
    return 4.0 * u[n][j] - u[n - 1][j];
}
void x_filling(double *x_points, int n, double start, double finish)
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
void file_write(FILE *file_graph, double **u, int time_parts, int space_parts, double *time_points, double *space_points)
{
    FILE *f1 = fopen("dataut.txt", "w");

    FILE *f2 = fopen("dataux.txt", "w");
    for (int n = 0; n < time_parts; ++n)
    {
        for (int j = 0; j < space_parts; ++j)
            fprintf(file_graph, "%.15lf %.15lf %.15lf\n", time_points[n], space_points[j], u[n][j]);
    }
    for (int i = 0; i < time_parts; ++i)
        fprintf(f1, "%.15lf %.15lf\n", time_points[i], u[i][(int)(space_parts-1)]);
    for (int i = 0; i < space_parts; ++i)
        fprintf(f2, "%.15lf %.15lf\n", space_points[i], u[(int)(time_parts-1)][i]);
}
// double max_error(double *points_x, double *points_y, int n)
// {
//     double max = 0;
//     for (int i = 0; i < n; ++i)
//     {
//         if (fabs(points_y[i] - utoch(points_x[i])) > max)
//             max = fabs(points_y[i] - utoch(points_x[i]));
//     }
//     return max;
// }
double find_max(int time_parts, int space_parts, double **u, double **u_touch)
{
    double max_value = 0;
    double sred = 0;
    for(int i = 0; i < time_parts; ++i)
    {
        for(int j = 0; j < space_parts; ++j)
        {
            if(fabs(u_touch[i][j] - u[i][j]) >= max_value)
            {
                max_value = fabs(u_touch[i][j] - u[i][j]);
            }
            sred+=u[i][j];
        }
    }
    printf("%.15lf\n", sred/(space_parts*time_parts));
    return max_value;
}
void progonka(double *time_points, double *space_points, double **u, int time_parts, int space_parts, double time_step, double space_step, int n)
{
    double *k_points = malloc(sizeof(double) * space_parts);
    double *v_points = malloc(sizeof(double) * space_parts);
    k_points[0] = 0.0;
    v_points[0] = 1.0;
    for (int i = 1; i < space_parts; ++i)
    {
        k_points[i] = -(Bi(time_step, space_step) / (Ci(time_step, space_step) + Ai(time_step, space_step) * k_points[i - 1]));
    }
    for (int i = 1; i < space_parts; ++i)
    {
        v_points[i] = -((Ai(time_step, space_step) * v_points[i - 1] - Fi(n - 1, i, u)) / (Ci(time_step, space_step) + Ai(time_step, space_step) * k_points[i - 1]));
    }
    u[n][space_parts - 1] = (v_points[space_parts - 1] + k_points[space_parts - 1] * v_points[space_parts - 2]) / (1 - k_points[space_parts - 1] * k_points[space_parts - 2]);
    for (int i = space_parts - 2; i >= 0; --i)
    {
        u[n][i] = k_points[i] * u[n][i + 1] + v_points[i];
    }
}
void shema(double *time_points, double *space_points, double **u, int time_parts, int space_parts, double time_step, double space_step, int n)
{
     for (int i = 1; i < space_parts; ++i)
    {
        u[n][i] = u[n-1][i] - (time_step / space_step) * (u[n-1][i] - u[n-1][i - 1]);
    }
}
int main()
{
    int space_parts = 1000, time_parts = 1000;
    double space_start = 0.0, space_end = 1.0, time_start = 0.0, time_end = 1.0, time_step = (time_end - time_start) / (time_parts - 1);
    double space_step = (space_end - space_start) / (space_parts - 1);
    double **u = malloc(sizeof(double *) * time_parts);
    double **u_touch = malloc(sizeof(double *) * time_parts);
    double *time_points = malloc(sizeof(double) * time_parts);
    double *space_points = malloc(sizeof(double) * space_parts);
    x_filling(time_points, time_parts, time_start, time_end);
    x_filling(space_points, space_parts, space_start, space_end);
    FILE *f = fopen("data.txt", "w");
    for (int i = 0; i < time_parts; ++i)
    {
        u[i] = malloc(sizeof(double) * space_parts);
    }
    for (int i = 0; i < time_parts; ++i)
    {
        u_touch[i] = malloc(sizeof(double) * space_parts);
    }
    for(int i = 0; i<time_parts;++i)
    {
        for(int j = 0; j < space_parts; ++j)
        {
            if(time_points[i] > space_points[j])
                u_touch[i][j] = 0;
            else
            {
                u_touch[i][j] = 1;
            }
            
        }
    }
    for (int i = 0; i < time_parts; ++i)
    {
        u[i][0] = 1;
    }
    for (int i = 0; i < space_parts; ++i)
    {
        u[0][i] = 0;
    }
    for (int i = 1; i < space_parts; ++i)
    {
        u[1][i] = u[0][i] - (time_step / space_step) * (u[0][i] - u[0][i - 1]);
    }
    for (int n = 2; n < time_parts; ++n)
    {
        shema(time_points, space_points, u, time_parts, space_parts, time_step, space_step, n);
    }
    file_write(f, u, time_parts, space_parts, time_points, space_points);
    printf("%.15lf\n", find_max(time_parts, space_parts, u, u_touch));
    return 0;
}
