#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "coeffs.h"
#include "coeffs_EZ.h"

using namespace std;

#define BUFF_SIZE 136
/*extern float rawECG[BUFF_SIZE];
extern float notchECG[BUFF_SIZE];
extern float filtECG[BUFF_SIZE];
extern float derECG[BUFF_SIZE];
extern float sqECG[BUFF_SIZE];
extern float movECG[BUFF_SIZE];
extern float lfiltECG[BUFF_SIZE];*/

extern float signal[BUFF_SIZE]; //store signal segment to be filtered
extern float bp_signal[BUFF_SIZE];
extern float derivative[BUFF_SIZE];
extern float squared[BUFF_SIZE];
extern float integral[BUFF_SIZE];
extern bool outputSignal[BUFF_SIZE];
extern float r_peaks[BUFF_SIZE];

/*Final heart rate stored here*/
extern float HR[BUFF_SIZE];

extern long unsigned int current;



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

void lfilter()
{

  for(int i=0; i < BUFF_SIZE; i++)
	{
		float tmp = 0.;
		int j=0;
		lfiltECG[i] = 0.f;
		for(j=0; j < 9; j++)
		{
			if(i - j < 0) continue;
			tmp += EZ_B[j] * movECG[i-j];
		}

		for(j=1; j < 9; j++)
		{
			if(i - j < 0) continue;
			tmp -= EZ_A[j]*lfiltECG[i-j];
		}

		tmp /= EZ_A[0];
		lfiltECG[i] = tmp;
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

void diffFunc(float *in)
{
  
  for(int i = 1; i < BUFF_SIZE;i++)
  {
    HR[i-1] = (in[i]-in[i-1]);
    HR[i-1] = HR[i-1]/200;
    HR[i-1] = 60/HR[i-1];
  }

 
}

void ez_detector_iter(double *sig, int sig_length, int fs)
{
  /*double f1 = 48/fs;
  double f2 = 52/fs;
  lfilter();

  int count = 0; 
  double diff[BUFF_SIZE];
  diff[0] = 0.0;
  diff[1] = 0.0;
  diff[2] = 0.0;
  diff[3] = 0.0;
  for(int i = 4; i < sig_length; i++)
  {
    diff[i] = lfiltECG[i] - lfiltECG[i-4];
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
  }*/


  double m200 = 0.2*fs;
  double m1200 = 1.2*fs;
  double m160 = 0.16*fs;
  int ms200 = (int)m200;
  int ms1200 = (int)m1200;
  int ms160 = (int)m160;

  double neg_thresh = 0.01*fs;
  int neg_threshold = (double)neg_thresh;

  //adaptive steep slope threshold
  double M = 0;
  
  //A buffer with 5 steep slope threshold values
  vector<double> MM;
  
  //M5 (an element of MM) is recaclculated periodically based on the sampling frequency
  double newM5 = 0;
  
  //
  vector<double> M_list;
  vector<double> neg_m;
  
  vector<double> MM;
  //need m_slope

  vector<double> QRS;

  int counter = 0;
  vector<double> thi_list;
  bool thi = false;
  vector<double> thf_list;
  bool thf = false;

  for(int i = 0; i < BUFF_SIZE; i++)
  {
    if(i < 5*fs)
    {
      M = 0.6*max_vec(integral, 0, i+1);
      MM.push_back(M);
      if(MM.size() > 5)
      {
        MM.erase(MM.begin());
      }
    }

      else if(QRS.size() > 0 && i < QRS[QRS.size()-1]+ms200)
      {
        newM5 = 0.6*max_vec(integral, QRS[QRS.size()-1], i);
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

      if(QRS.size() == 0  && integral[i] > M)
      {
        QRS.push_back(i);
        thi_list.push_back(i);
        thi = true;
      }

      else if(QRS.size()>0 && i > QRS[QRS.size()-1]+ms200 && integral[i]>M)
      {
        QRS.push_back(i);
        thi_list.push_back(i);
        thi = true;
      }

      if(thi && i <thi_list[thi_list.size()-1]+ms160)
      {
        if(integral[i]<-M && integral[i-1] > -M)
        {
            thf = true;
        }
        if(thf && integral[i] < -M)
        {
          thf_list.push_back(i);
          counter++;
        }
        else if(integral[i] > -M && thf)
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
      r_peaks[count] = push_val;
	count++;
      counter = 0;
      thi = false;
      thf = false;
      if(count == BUFF_SIZE-1) break;
    }


  }
  
   diffFunc(r_peaks);

}

