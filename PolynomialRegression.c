#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void LUdecompose(int n, double iarr[n][n], double Uarr[n][n], double Larr[n][n])
{	
	//initialize to 0
	for(int y = 0; y < n; y++)
	{
		for(int x = 0; x < n; x++)
		{
			Larr[x][y] = 0;
			Uarr[x][y] = 0;
		}
	}

	for(int x = 0; x < n; x++)
	{
		for(int z = x; z < n; z++) //up
		{
			double sum = 0;
			for(int y = 0; y < x; y++)
			{
				sum += Larr[y][x] * Uarr[z][y];
			}
			Uarr[z][x] = iarr[z][x] - sum;
		}
		for(int z = x; z < n; z++) //low
		{
			if(x == z) //diag
			{
				Larr[x][x] = 1;
				continue;
			}
			double sum = 0;
			for(int y = 0; y < x; y++)
			{
				sum += Larr[y][z] * Uarr[x][y];
			}
			Larr[x][z] = (iarr[x][z] - sum)/Uarr[x][x];
		}
	}
}
void SolveLU(int n, double Xarr[n], double Uarr[n][n], double Larr[n][n], double Rarr[n])
{
    //LUx = b; Ux = y
	//solve y for Ly = b
	double Yarr[n];
	for(int y = 0; y < n; y++)//initialize
	{
		Yarr[y] = 0;
	}
	for(int y = 0; y < n; y++)
	{
		double sum = 0;
		for(int x = 0; x < y; x++)
		{	
			(x != y)?(sum += Larr[x][y] * Yarr[x]):printf("");
		}
		Yarr[y] = (Rarr[y] - sum)/Larr[y][y];
	}
	//LUx = b; Ux = y
	//solve x for Ux = y
	for(int y = 0; y < n; y++)//initialize
	{
		Xarr[y] = 0;
	}
	for(int y = n-1; y >= 0; y--)
	{
		double sum = 0;
		for(int x = n-1; x >= y; x--)
		{
			(x != y)?(sum += Uarr[x][y] * Xarr[x]):printf("");
		}
		(Uarr[y][y] == 0)?(Xarr[y] = 0):(Xarr[y] = (Yarr[y] - sum)/Uarr[y][y]);
	}
}
double SigmaFunc1(int n, double xarr[n], int power)
{
    double result = 0;
    for(int x = 0; x < n; x++)
    {
        result += pow(xarr[x], power);
    }

    return result;
}
double SigmaFunc2(int n, double xarr[n], double yarr[n], double power)
{
    double result = 0;
    for(int x = 0; x < n; x++)
    {
        result += pow(xarr[x], power)*yarr[x];
    }

    return result;
}
void PolynomialRegression(int n, double xarr[n], double yarr[n], int r, double Parr[r])
{
    //generate matrix
    double Carr[r][r]; //coefficient
    for(int x = 0; x < r; x++)
    {
        for(int y = 0; y < r; y++)
        {
            if(x == 0 && y == 0)
            {
                Carr[x][y] = n;
            }
            else
            {
                Carr[x][y] = SigmaFunc1(n, xarr, x+y);
            }
        }
    }
    double Rarr[r]; //result
    for(int x = 0; x < r; x++)
    {
        Rarr[x] = SigmaFunc2(n, xarr, yarr, x);
    }

    double Uarr[r][r]; //upper
    double Larr[r][r]; //lower
    LUdecompose(r, Carr, Uarr, Larr);
    SolveLU(r, Parr, Uarr, Larr, Rarr);
}
double CalcPolynomial(double x_val, int r, double coefficients[r])
{
    double result = 0;
    for(int x = 0; x < r; x++)
    {
        result += coefficients[x]*pow(x_val, x); 
    }
    return result;
}
double CorrelationCoefficient(int n, double xarr[n], double yarr[n], int r, double coefficients[r])
{
    //y_avg
    double y_avg;
    double y_total = 0;
    for(int x = 0; x < n; x++)
    {
        y_total += yarr[x];
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
        D += pow((yarr[x]-(CalcPolynomial(xarr[x], r, coefficients))),2);
    }
	printf("(Dt=%lf;D=%lf)\n", Dt, D);
    return(sqrt((Dt-D)/Dt));
}
void calc(int data_n, double Xarr[data_n], double Yarr[data_n], int r)
{
	printf("n=%d\n", data_n);
    double Carr[r];
    PolynomialRegression(data_n, Xarr, Yarr, r, Carr); 

	printf("fungsi regresi polynomial:\n", r-1);
	for(int y = 0; y < r; y++)
	{
		if(y!=0 && Carr[y]>0)
		{
			printf("+");
		}
		printf("%lf*(t^%d)", Carr[y], y);
	}
	printf("\n");
	double correlation_coefficient = CorrelationCoefficient(data_n, Xarr, Yarr, r, Carr);
    printf("r = %lf", correlation_coefficient);
	printf("\nr^2 = %lf", pow(correlation_coefficient,2));
}
void main()
{
	int r = 4;
    double xAxis[128];
    for(int x = 0; x < 128; x++)
    {
        xAxis[x] = x+1;
    }

    double Data[] = {554,771,1208,1870,2613,4349,5739,7417,9308,11289,13748,16369,19383,22942,26302,28985,31774,33738,35982,37626,38791,51591,55748,56873,57416,57934,58016};
    calc((sizeof(Data)/sizeof(double)), xAxis, Data, r);
}
