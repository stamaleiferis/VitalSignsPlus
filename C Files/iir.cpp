#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "coeffs.h"
#include "coeffs_EZ.h"

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

void lfilter(double *b, double *a, double *in, double *out, int in_length, int a_length, int b_length)
{
  for(int i=0; i < in_length; i++)
	{
		float tmp = 0.;
		int j=0;
		out[i] = 0.f;
		for(j=0; j < b_length; j++)
		{
			if(i - j < 0) continue;
			tmp += b[j] * in[i-j];
		}

		for(j=1; j < a_length; j++)
		{
			if(i - j < 0) continue;
			tmp -= a[j]*out[i-j];
		}

		tmp /= a[0];
		out[i] = tmp;
	}
}
double mean(vector<double> MM)
{
double mean = 0;
for(int i = 0; i < MM.size(); i++)
{
  mean += MM[i];
}
mean/=MM.size();
return mean;
}

double max_vec(vector<double> arr, int start, int end)
{
  double max = 0;
  for(int i = start; i < end; i++)
  {
    if(arr[i]> max)
    {
      max = arr[i];
    }
  }
  return max;
}

int argmax(vector<double> arr, int start, int end)
{
  int max = 0;
  for(int i = start; i < end; i++)
  {
    if(arr[i] > arr[max])
    {
      max = i;
    }
  }
  return max;
}
vector<int> engzee_detector(double *sig, int sig_length, int fs)
{
  double f1 = 48/fs;
  double f2 = 52/fs;
  double output[sig_length];
  lfilter(EZ_B, EZ_A, sig, output, sig_length, 9, 9);

  double diff[sig_length];
  diff[0] = 0.0;
  diff[1] = 0.0;
  diff[2] = 0.0;
  diff[3] = 0.0;
  for(int i = 4; i < sig_length; i++)
  {
    diff[i] = output[i] - output[i-4];
  }



  double c1[] = {1.0, 4.0, 6.0, 4.0, 1.0};
  double scalar_1[] = {1.0};
  double lp[sig_length+4];
  conv(diff, c1, lp, sig_length, 5);
  vector<double> low_pass;
  for(int i = 0; i < sig_length; i++)
  {
    low_pass.push_back(lp[i]);
  }
  int stop = 0.2*fs;
  for(int i = 0; i < stop; i++)
  {
    low_pass[i] = 0;
  }


  double m200 = 0.2*fs;
  double m1200 = 1.2*fs;
  double m160 = 0.16*fs;
  int ms200 = (int)m200;
  int ms1200 = (int)m1200;
  int ms160 = (int)m160;

  double neg_thresh = 0.01*fs;
  int neg_threshold = (double)neg_thresh;

  double M = 0;
  double newM5 = 0;
  vector<double> M_list;
  vector<double> neg_m;
  vector<double> MM;
  //need m_slope

  vector<double> QRS;
  vector<int> r_peaks;

  int counter = 0;
  vector<double> thi_list;
  bool thi = false;
  vector<double> thf_list;
  bool thf = false;

  for(int i = 0; i < low_pass.size(); i++)
  {
    if(i < 5*fs)
    {
      M = 0.6*max_vec(low_pass, 0, i+1);
      printf("%f\n", max_vec(low_pass, 0, i+1));
      MM.push_back(M);
      if(MM.size() > 5)
      {
        MM.erase(MM.begin());
      }
    }

      else if(QRS.size() > 0 && i < QRS[QRS.size()-1]+ms200)
      {
        newM5 = 0.6*max_vec(low_pass, QRS[QRS.size()-1], i);
        if(newM5 > 1.5*MM[MM.size()-1])
        {
          newM5 = 1.1*MM[MM.size()-1];
        }
      }

      else if(QRS.size() > 0 && i == QRS[QRS.size()-1]+ms200)
      {
        MM.push_back(newM5);
        if(MM.size() > 5)
        {
          MM.erase(MM.begin());
        }

        M = mean(MM);
      }

      else if(QRS.size() > 0 && i > QRS[QRS.size()-1]+ms200 && i < QRS[QRS.size()-1]+ms1200)
      {
        M = mean(MM);
        int ind = QRS[QRS.size()-1]+ms200;
        M*=M_slope[i-ind];
      }

      else if(QRS.size()>0 && i> QRS[QRS.size()-1]+ms1200)
      {
        M = 0.6*mean(MM);
      }

      M_list.push_back(M);
      neg_m.push_back(-M);

      if(QRS.size() == 0  && low_pass[i] > M)
      {
        QRS.push_back(i);
        thi_list.push_back(i);
        thi = true;
      }

      else if(QRS.size()>0 && i > QRS[QRS.size()-1]+ms200 && low_pass[i]>M)
      {
        QRS.push_back(i);
        thi_list.push_back(i);
        thi = true;
      }

      if(thi && i <thi_list[thi_list.size()-1]+ms160)
      {
        if(low_pass[i]<-M && low_pass[i-1] > -M)
        {
            thf = true;
        }
        if(thf && low_pass[i] < -M)
        {
          thf_list.push_back(i);
          counter++;
        }
        else if(low_pass[i] > -M && thf)
        {
          counter = 0;
          thi = false;
          thf = false;
        }
      }
      else if(thi && i > thi_list[thi_list.size()-1]+ms160)
      {
        counter = 0;
        thi = false;
        thf = false;
      }

      if(counter > neg_threshold)
      {
      int start = thi_list[thi_list.size()-1]-(int)0.01*fs;
      vector<double> unfiltered_section;
      for(int j = start; i < i; j++)
      {
        unfiltered_section.push_back(sig[j]);
      }

      double x = 0.01*fs;
      int x_push = (int)x;
      int push_val = argmax(unfiltered_section, 0, unfiltered_section.size())+thi_list[thi_list.size()-1]-x;
      r_peaks.push_back(push_val);
      counter = 0;
      thi = false;
      thf = false;
    }


  }
  return r_peaks;


}

vector<double> diff(vector<int> in)
{
  vector<double> differences;
  for(int i = 1; i < in.size();i++)
  {
    differences.push_back(in[i]-in[i-1]);
    differences[i-1] = differences[i-1]/200;
    differences[i-1] = 60/differences[i-1];
  }
  return differences;
}

int main(double *raw_sig)
{
    double sig[i-300];
    for(int j = 300; j < i; j++)
    {
      sig[j-300] = raw_sig[j];
    }

    /*60 Hz, 60 Hz, 80 Hz notch filtering*/
    int length = i;
    double out[length];
    double out_60[length];
    double out_80[length];
    int filter_length = 3;
    cur_length = length;
    filter1(NOTCH_60_B, NOTCH_60_A, sig, out, length, filter_length);
    filter1(NOTCH_60_B, NOTCH_60_A, out, out_60, length, filter_length);
    filter1(NOTCH_80_B, NOTCH_80_A, out_60, out_80, length, filter_length);

    /*fir filtering*/
    double fir[i+128-1];
    double fir2[i+128-1];
    conv(out_80, TAPS_KAISER10, fir, i, 128);
    conv(fir, TAPS_KAISER16,fir2, i, 128);

    /*derivative filter*/
    double der_h[5] = {-0.125, -0.25, 0, 0.25, 0.125};
    double derecg[length+4];
    conv(fir2, der_h, derecg, length, 5);
    float max = max_index(derecg, length+4);
    double sqECG[length+4];
    for(int j = 0; j < length+4; j++)
    {
      sqECG[j] = derecg[j] / max;
      sqECG[j] = pow(sqECG[j], 2.0);
    }

    double h[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    double movECG[length+13];
    conv(sqECG, h, movECG, length, 10);

    /*call to EZ, returns peaks indices*/
    vector<int> peaks = engzee_detector(movECG, length+13, 200);

    /*uses peak indices to get heart rates*/
    vector<double> HR = diff(peaks);


}
