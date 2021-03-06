#ifndef RECORDINDICES_H
#define RECORDINDICES_H

#define iVar 1
#define iLower 2
#define iUpper 3
#define iStep 4
#define iAdmin 1
#define FHRecordR 2
#define iinvAlfaMZ 2
#define iAlfasMZ 3
#define iGF 4
#define iME 5
#define iMU 6
#define iMD 7
#define iMM 8
#define iMC 9
#define iMS 10
#define iML 11
#define iMT 12
#define iMB 13
#define iMW 14
#define iMZ 15
#define iCKMlambda 16
#define iCKMA 17
#define iCKMrhobar 18
#define iCKMetabar 19
#define iTB 20
#define iMA0 21
#define iMHp 22
#define iMSusy 23
#define iM1SL 24
#define iM1SE 25
#define iM1SQ 26
#define iM1SU 27
#define iM1SD 28
#define iM2SL 29
#define iM2SE 30
#define iM2SQ 31
#define iM2SU 32
#define iM2SD 33
#define iM3SL 34
#define iM3SE 35
#define iM3SQ 36
#define iM3SU 37
#define iM3SD 38
#define iQtau 39
#define iQt 40
#define iQb 41
#define iscalefactor 42
#define iprodSqrts 43
#define FHRecordC 44
#define iAe 44
#define iAu 48
#define iAd 52
#define iAmu 56
#define iAc 60
#define iAs 64
#define iAtau 68
#define iAt 72
#define iAb 76
#define iXtau 80
#define iXt 84
#define iXb 88
#define iMUE 92
#define iM1 96
#define iM2 100
#define iM3 104
#define ideltaLLL12 108
#define ideltaLLL23 112
#define ideltaLLL13 116
#define ideltaELR12 120
#define ideltaELR23 124
#define ideltaELR13 128
#define ideltaERL12 132
#define ideltaERL23 136
#define ideltaERL13 140
#define ideltaERR12 144
#define ideltaERR23 148
#define ideltaERR13 152
#define ideltaQLL12 156
#define ideltaQLL23 160
#define ideltaQLL13 164
#define ideltaULR12 168
#define ideltaULR23 172
#define ideltaULR13 176
#define ideltaURL12 180
#define ideltaURL23 184
#define ideltaURL13 188
#define ideltaURR12 192
#define ideltaURR23 196
#define ideltaURR13 200
#define ideltaDLR12 204
#define ideltaDLR23 208
#define ideltaDLR13 212
#define ideltaDRL12 216
#define ideltaDRL23 220
#define ideltaDRL13 224
#define ideltaDRR12 228
#define ideltaDRR23 232
#define ideltaDRR13 236
#define FHRecordE 240
#define FHRecordN 239

#endif
* FHRecord.h.in
* the data structures for a FH record
* this file is part of FeynHiggs
* last modified 11 Jul 11 th


#ifndef FHRangeR
#define FHRangeR FHRecordR, FHRecordC - 1
#define FHRangeC FHRecordC, FHRecordN, 4
#define FHRangeA FHRecordR, FHRecordN

#define iRe(i) i
#define iIm(i) i+1
#define iAbs(i) i+2
#define iArg(i) i+3

#define iMSS(n,g) iM1SL+(n-1)*(iM1SE-iM1SL)+(g-1)*(iM2SL-iM1SL)
#define iMf(t,g) iMU+(t-3)*(iMD-iMU)+(g-1)*(iMC-iMU)
#define iAf(t,g) iAu+(t-3)*(iAd-iAu)+(g-1)*(iAc-iAu)
#define iQSf(t) iQt+(t-3)*(iQb-iQt)

#define FHNameR(i) FHName(i)(1:len_trim(FHName(i)))
#define FHNameC(i) FHName(i)(4:index(FHName(i),")")-1)

#define RecordDecl(rec) double precision rec(FHRecordN,4)
#endif

	double precision unset, default, bytable
	parameter (unset = -999)
	parameter (default = -888)
	parameter (bytable = 777)

	character*16 FHName(FHRecordR:FHRecordN)
	common /fhrecnames/ FHName

	integer maxcols, maxrows
	parameter (maxcols = FHRecordN, maxrows = 2400)

	double precision tabledata(maxcols,maxrows)
	integer tableflag(0:maxcols), tablerows
	common /fhtable/ tabledata, tableflag, tablerows

