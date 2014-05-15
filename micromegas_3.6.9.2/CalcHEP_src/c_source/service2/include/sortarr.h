#ifndef SORTARR
#define SORTARR(arr,len) {int Ind_=1;long L;  while (Ind_ < (len)) \
if (arr[Ind_-1] <= arr[Ind_]) Ind_++; else {L=arr[Ind_-1];arr[Ind_-1]=arr[Ind_];arr[Ind_]=L;if(Ind_>1)Ind_--;}} 
#endif
