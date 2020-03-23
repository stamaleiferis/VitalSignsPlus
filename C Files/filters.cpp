#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "coeffs.h"
#include "coeffs_EZ.h"

using namespace std;

extern float rawECG[BUFF_SIZE];
extern float notchECG[BUFF_SIZE];

extern float filtECG[BUFF_SIZE];
extern float derECG[BUFF_SIZE];
extern float sqECG[BUFF_SIZE];
extern float movECG[BUFF_SIZE];

void filter_helper(double *coeffs_B, double *coeffs_A, double *input, double *output, int length, int filterLength)
{
    double bcc, acc;
    double *inputp;
    int n,k;
    for (int ii=0; ii<filterLength; ii++)
    {
        output[ii] = 0;
    }
    output[0] = coeffs_B[0] * input[0];
    output[1] = (coeffs_B[0] * input[1] + coeffs_B[1] * input[0]) - coeffs_A[1] * output[0];

    for (n = 2; n < length; n++) {

        acc = 0;
        bcc = 0;

        for (k = 0; k < filterLength; k++)
        {
            bcc += coeffs_B[k] * input[n-k]; //b[0] * x[6] + b[1] * x[5]...+b[6] * x[0]
        }

        for (k = 0; k < filterLength-1; k++)
        {
            acc += coeffs_A[k] * output[n-k-1]; //a[1] * y[5] + a[2] * y[4]...+a[6] * y[0]
        }
        output[n] = bcc-acc;

    }


}

void notchFilters()
{

  double notchECG_1[BUFF_SIZE];
  double notchECG_2[BUFF_SIZE];
  //apply 60Hz Notch Filter
  filter_helper(NOTCH_60_B, NOTCH_60_A, rawECG, notchECG_1, BUFF_SIZE, 3);

  //apply 60Hz Notch Filter
  filter_helper(NOTCH_60_B, NOTCH_60_A, notchECG_1, notchECG_2, BUFF_SIZE, 3);

  //apply 80Hz Notch Filter
  filter_helper(NOTCH_60_B, NOTCH_60_A, notchECG_2, notchECG, BUFF_SIZE, 3);

}

void conv(double *A, double *B, double *out, int lenA, int lenB)
{
	int nconv;
	int i, j, i1;
	float tmp;
	//allocated convolution array
	nconv = lenA+lenB-1;
	//convolution process
  int count = 0;
	for (i=0; i<nconv; i++)
	{
		i1 = i;
		tmp = 0.0;
		for (j=0; j<lenB; j++)
		{
			if(i1>=0 && i1<lenA)
				tmp = tmp + (A[i1]*B[j]);
			i1 = i1-1;
		}
    out[count++] = tmp;
	}

}

void kaiserFilters()
{
  double fir1[BUFF_SIZE];
  conv(notchECG, TAPS_KAISER10, fir1, BUFF_SIZE, 128);
  conv(fir1, TAPS_KAISER16, filtECG, BUFF_SIZE, 128);
}

void derivativeFilters()
{
  double der_h[5] = {-0.125, -0.25, 0, 0.25, 0.125};
  conv(filtECG, der_h, derECG, BUFF_SIZE, 5);
}

float max(double *arr, int arr_length)
{
  int max = 0;
  for(int i = 0; i < arr_length; i++)
  {
    if(arr[i] > arr[max])
    {
      max = i;
    }
  }
  return arr[max];
}

void squaring()
{
  float max = max(derECG, BUFF_SIZE);
  for(int j = 0; j < BUFF_SIZE; j++)
  {
    sqECG[j] = derECG[j] / max;
    sqECG[j] = pow(sqECG[j], 2.0);
  }

}

void movFilter()
{
  double h[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  conv(sqECG, h, movECG, BUFF_SIZE, 10);
}

int filtIter()
{
    /*60 Hz, 60 Hz, 80 Hz notch filtering
      output stored in notchECG*/
    notchFilters();

    /*fir filtering
      output stored in filtECG*/
    kaiserFilters();

    /*derivative filter
      output stored in derECG*/
    derivativeFilters();

    /*output stored in sqECG*/
    squaring();

    /*output stored in movECG*/
    movFilter();
}
