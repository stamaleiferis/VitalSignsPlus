#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "coeffs.h"

using namespace std;

int cur_length = 0;

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
    bool flipped = false;
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
void filter1(double *coeffs_B, double *coeffs_A, double *input, double *output, int length, int filterLength)
{
  filter_helper(coeffs_B, coeffs_A, input, output, length, filterLength);

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


float max_index(double *arr, int arr_length)
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

int main(int argc, char* argv[])
{
    if (argc != 2) perror("Incorrect num args");

    FILE *fp;
    fp = fopen(argv[1], "r");
    char buff[1024];
    int i = 0;
    while(fgets(buff, 1024, fp))
    {
      i++;
    }
    fclose(fp);
    fp = fopen(argv[1], "r");
    double sig[i];
    int count = 0;
    while(fgets(buff, 1024, fp))
    {
      sig[count++] = strtod(buff, NULL);
    }
    fclose(fp);


    int length = i;
    double out[length];
    double out_60[length];
    double out_80[length];
    int filter_length = 3;
    cur_length = length;
    filter1(NOTCH_60_B, NOTCH_60_A, sig, out, length, filter_length);
    filter1(NOTCH_60_B, NOTCH_60_A, out, out_60, length, filter_length);
    filter1(NOTCH_80_B, NOTCH_80_A, out_60, out_80, length, filter_length);

    double fir[i+128-1];
    double fir2[i+128-1];
    conv(out_80, TAPS_KAISER10, fir, i, 128);
    conv(fir, TAPS_KAISER16,fir2, i, 128);

    double der_h[5] = {-0.125, -0.25, 0, 0.25, 0.125};
    double derecg[length+4];
    conv(fir2, der_h, derecg, length, 5);
    float max = max_index(derecg, length+4);
    printf("%f\n", max);
    double sqECG[length+4];
    for(int j = 0; j < length+4; j++)
    {
      sqECG[j] = derecg[j] / max;
      sqECG[j] = pow(sqECG[j], 2.0);
    }

    double h[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double movECG[length+13];
    conv(sqECG, h, movECG, length, 10);


   fp = fopen("t.txt", "w");
   for(int j = 0; j < length+13; j++)
   {
     char buff[32];
     sprintf(buff, "%f\n", movECG[j]);
     fputs(buff, fp);
   }
   fclose(fp);

}
