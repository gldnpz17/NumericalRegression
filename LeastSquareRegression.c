#include<stdio.h>
#include<math.h>

double SigmaFunc1(int n, double xarr[n], double yarr[n])
{
    double result = 0;
    for(int x = 0; x < n; x++)
    {
        result += xarr[x]*yarr[x];
    }
    return result;
}
double SigmaFunc2(int n, double iarr[n])
{
    double result = 0;
    for(int x = 0; x < n; x++)
    {
        result += iarr[x];
    }
    return result;
}
double SigmaFunc3(int n, double xarr[n])
{
    double result = 0;
    for(int x = 0; x < n; x++)
    {
        result += xarr[x]*xarr[x];
    }
    return result;
}
double CalcCoefficient(int n, double xarr[n], double yarr[n])
{
    double result = (n*SigmaFunc1(n, xarr, yarr)-SigmaFunc2(n, xarr)*SigmaFunc2(n, yarr))/(n*SigmaFunc3(n, xarr)-SigmaFunc2(n, xarr)*SigmaFunc2(n, xarr));
    return result;
}
double CalcConstant(int n, double xarr[n], double yarr[n], double coefficient)
{
    double x_total = 0;
    double y_total = 0;
    for(int x = 0; x < n; x++)
    {
        x_total += xarr[x];
        y_total += yarr[x];
    }
    return (y_total/n - coefficient*(x_total/n));
}
double CorrelationCoefficient(int n, double xarr[n], double yarr[n], double coefficient, double constant)
{
    //y_avg
    double y_avg;
    double y_total = 0;
    for(int x = 0; x < n; x++)
    {
        y_total+=yarr[x];
    }
    y_avg = y_total/n;
    //Dt
    double Dt = 0;
    for(int x = 0; x < n; x++)
    {
        Dt += pow((yarr[x]-y_avg),2);
    }
    //D
    double D = 0;
    for(int x = 0; x < n; x++)
    {
        D += pow((yarr[x]-(constant+xarr[x]*coefficient)),2);
    }

    return(sqrt((Dt-D)/Dt));
}
void main()
{
    int data_n = 128;
    double xAxis[data_n];
    for(int x = 0; x < data_n; x++)
    {
        xAxis[x] = x+1;
    }
    double Data[] = {554,771,1208,1870,2613,4349,5739,7417,9308,11289,13748,16369,19383,22942,26302,28985,31774,33738,35982,37626,38791,51591,55748,56873,57416,57934,58016};

    double coefficient = CalcCoefficient(data_n, xAxis, Data);
    printf("koefisien = %lf\n", coefficient);
    double constant = CalcConstant(data_n, xAxis, Data, CalcCoefficient(data_n, xAxis, Data));
    printf("konstanta = %lf\n", constant);
    printf("r = %lf", CorrelationCoefficient(data_n, xAxis, Data, coefficient, constant));
}
