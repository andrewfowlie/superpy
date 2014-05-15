#include "stdafx.h"

#include <math.h>

#include "recipies.h"


float chebev(float a, float b, float c[], int m, float x)
{
	float d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) nrerror("x not in range in routine CHEBEV");
	y2=(float) 2.0*(y=((float) 2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return (float) (y*d-dd+0.5*c[0]);
}
