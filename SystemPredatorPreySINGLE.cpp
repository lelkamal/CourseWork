#include <mpi.h>
#include <iostream>
#include <math.h>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <time.h>

using namespace std;

int n;
float h, eps;
float maX = -1.0f;
float max1 = -1.0f;
float max2 = -1.0f;
float miN = 20000.0f;
float min1 = 20000.0f;
float min2 = 20000.0f;

float x1[10000];
float x2[10000];

float y11[10000];
float y2[10000];

float z1[10000];
float z2[10000];
float a, b, c, d;

void Search_h1(float h, int i)
{

	float q1, q2, q3, q4, p1, p2, p3, p4;

	if (z1[i] > maX)maX = z1[i];
	if (x1[i] > max1)max1 = x1[i];
	if (y11[i] > max2)max2 = y11[i];

	if (z1[i] < miN)miN = z1[i];
	if (x1[i] < min1)min1 = x1[i];
	if (y11[i] < min2)min2 = y11[i];

	p1 = h * (a - c * y11[i]) * x1[i];
	q1 = h * (-b + d * x1[i]) * y11[i];
	p2 = h * (a - c * (y11[i] + q1 / 2)) * (x1[i] + p1 / 2);
	q2 = h * (-b + d * (x1[i] + p1 / 2)) * (x1[i] + q1 / 2);
	p3 = h * (a - c * (y11[i] + q2 / 2)) * (x1[i] + p2 / 2);
	q3 = h * (-b + d * (x1[i] + p2 / 2)) * (x1[i] + q2 / 2);
	p4 = h * (a - c * (y11[i] + q3)) * (x1[i] + p3);
	q4 = h * (-b + d * (x1[i] + p3)) * (x1[i] + q3);
	x1[i + 1] = x1[i] + (p1 + 2 * p2 + 2 * p3 + p4) / 6;
	y11[i + 1] = y11[i] + (q1 + 2 * q2 + 2 * q3 + q4) / 6;
	z1[i + 1] = z1[i] + h;

	if (z1[i + 1] > maX)maX = z1[i + 1];
	if (x1[i + 1] > max1)max1 = x1[i + 1];
	if (y11[i + 1] > max2)max2 = y11[i + 1];

	if (z1[i + 1] < miN)miN = z1[i + 1];
	if (x1[i + 1] < min1)min1 = x1[i + 1];
	if (y11[i + 1] < min2)min2 = y11[i + 1];

	
}

void Search_h2(float h, int i){

	float q1, q2, q3, q4, p1, p2, p3, p4;

	p1 = h * (a - c * y2[i]) * x2[i];
	q1 = h * (-b + d * x2[i]) * y2[i];
	p2 = h * (a - c * (y2[i] + q1 / 2)) * (x2[i] + p1 / 2);
	q2 = h * (-b + d * (x2[i] + p1 / 2)) * (x2[i] + q1 / 2);
	p3 = h * (a - c * (y2[i] + q2 / 2)) * (x2[i] + p2 / 2);
	q3 = h * (-b + d * (x2[i] + p2 / 2)) * (x2[i] + q2 / 2);
	p4 = h * (a - c * (y2[i] + q3)) * (x2[i] + p3);
	q4 = h * (-b + d * (x2[i] + p3)) * (x2[i] + q3);
	x2[i + 1] = x2[i] + (p1 + 2 * p2 + 2 * p3 + p4) / 6;
	y2[i + 1] = y2[i] + (q1 + 2 * q2 + 2 * q3 + q4) / 6;
	z2[i + 1] = z2[i] + h;

	p1 = h * (a - c * y2[i]) * x2[i];
	q1 = h * (-b + d * x2[i]) * y2[i];
	p2 = h * (a - c * (y2[i] + q1 / 2)) * (x2[i] + p1 / 2);
	q2 = h * (-b + d * (x2[i] + p1 / 2)) * (x2[i] + q1 / 2);
	p3 = h * (a - c * (y2[i] + q2 / 2)) * (x2[i] + p2 / 2);
	q3 = h * (-b + d * (x2[i] + p2 / 2)) * (x2[i] + q2 / 2);
	p4 = h * (a - c * (y2[i] + q3)) * (x2[i] + p3);
	q4 = h * (-b + d * (x2[i] + p3)) * (x2[i] + q3);
	x2[i + 1] = x2[i] + (p1 + 2 * p2 + 2 * p3 + p4) / 6;
	y2[i + 1] = y2[i] + (q1 + 2 * q2 + 2 * q3 + q4) / 6;
}


void Runge()
{
	int i = 0;
	int k = 0;
	while (i < n){
		Search_h1(h, i);
		Search_h2(h / 2, i);

		if (fabs(x1[i] - x2[i]) >= eps || fabs(y11[i] - y2[i]) >= eps)
		{
			h = h / 2;

			i = 0;
			max1 = -1.0f;
			max2 = -1.0f;
		}
		else i++;
	}
}


int main(int argc, char* argv[])
{
	double t1, t2, dt;

	int size,rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	t1 = MPI_Wtime();

    for(int v=0;v<4;++v)
    {
        for (int i = 0; i < n; i++)
        {
            x1[i] = 0;
            x2[i] = 0;
            y11[i] = 0;
            y2[i] = 0;
        }
    
        char buf[100];
        sprintf(buf, "input%d.lm", v);
        freopen(buf, "r", stdin);

        cin >> n;
        cin >> a >> b >> c >> d;
        cin >> x1[0] >> y11[0];

        fclose(stdin);

        x2[0] = x1[0];
        y2[0] = y11[0];

        x2[0] = x1[0];
        y2[0] = y11[0];

        h = 0.1;
        eps = 0.001;

        Runge();

	
        sprintf(buf, "Soutput%d.lm", v);
        freopen(buf, "w", stdout);

        for (int j = 0; j < n; j++)
        {
            cout << z1[j] << " " << x1[j] << " " << y11[j] << endl;
        }

        fclose(stdout);
    }

	t2 = MPI_Wtime();

	dt = t2 - t1;

	freopen("TimeS.lm", "w", stdout);
	cout << dt << endl;
	fclose(stdout);

	MPI_Finalize();
	
	return 0;
}


