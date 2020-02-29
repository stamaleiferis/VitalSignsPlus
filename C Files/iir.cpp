#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;


void filter1(double *coeffs_B, double *coeffs_A, double *input, double *output, int length, int filterLength)
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
    FILE *fp = fopen("t.txt", "w");
    char buff[1024];
    for(int i = 0; i < n; i++)
    {
          sprintf(buff, "%f\n", output[i]);
          fputs(buff, fp);
    }
    fclose(fp);
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
    double out[i];
    double b[3] = {0.97697628, 0.60380455, 0.97697628};
    double a[3] = {0.60380455, 0.95395256};
    int length = i;
    int filter_length = 3;

    filter1(b, a, sig, out, length, filter_length);


}
