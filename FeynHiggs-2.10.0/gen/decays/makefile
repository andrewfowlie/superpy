all: sm mssm

ferm: smferm mssmferm
bos: smbos mssmbos

sm: smbos smferm
mssm: mssmbos mssmferm

smbos: f/hgagaSM f/hgaZSM f/hggSM
smferm: f/heeSM f/huuSM f/hddSM

mssmbos: f/hgaga f/hgaZ f/hgg f/hhh
mssmferm: f/hee f/huu f/hdd

f/%:
	./decay$(SUFFIX) $@

