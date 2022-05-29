#!/usr/bin/python3
#adapted by 'QUMU00'
import sys

finput = sys.argv[1]
foutput = sys.argv[2]
#NT='TCGA'.split('')
NT = ['T','C','G','A']
kmer = ['T','C','G','A']
mer = ['T','C','G','A']
NTT = ['T','C','G','A','[TCGA]']
#print (NTT)
NTN = ['T','C','G','A','[TCGA]']

tmp=[]
for AA in NT:
  add=AA
  for BB in mer:
    strr = BB+add
    tmp.append(strr)
mer = tmp
kmer = kmer + tmp  

#print(NT)
#for (i=3;i<6;i=i+1):
for i in range(3,6,1): 
  for AA in NT:
    add = AA
    for BB in NTN:	
      strr = ''
      strr = add+BB
      #print(strr)
      for CC in NT:
        strr2 = '' 
        #print(strr2,CC)
        strr2 = strr+CC
        kmer.append(strr2)
  
  tmp=[]
  for AA in NTN:
    strr = AA
    for BB in NTT:
      tmpp = ''
      tmpp = strr + BB
      tmp.append(tmpp)
  
  NTN=tmp

def mmm (f,ss):
  ff=f.replace('[TCGA]','N')
  for ii in range(0,len(ff)):
    if(ff[ii]=='N'):
      continue
    if(ff[ii]!=ss[ii]):
      return 0

  return 1 
    

#print (kmer)
ID=open(finput,'r+')
OUT=open(foutput,'w+')
strr=''
for ID_line in ID.readlines():
  strr = ID_line
  strr = strr.strip()
  strr = strr.upper()
  n = 0
  #print(kmer)
  for AA in kmer:
    #print(AA)   
    if('[' not in AA):
      n=n+1
      num = strr.count(AA)
      #if(!num):
      #  num=0
      lenn = len(AA)
      feq = num / (len(strr)-lenn + 1)
      OUT.write('%f\t'%feq)
      #print('%s,%d,%s,%f\t'%(strr,num,AA,feq))
    else:
      f=''
      lenn = 0
      f = AA
      n = n+1
      if(n > 20 and n<101):
        lenn = 3
      elif (n < 501):
        lenn = 4
      else:
        lenn = 5
      end = len(strr) - lenn
      num = 0
      #for ( i= 0: i<= end; i=i+1):
      for i in range(0,end+1,1):
        ss = strr[i:i+lenn]
        #if(f in ss ):
        if(mmm(f,ss)):
          num=num+1
      feq = num / (len(strr) - lenn + 1)
      OUT.write('%f\t'%feq)
      #print('%f\t'%feq)
  OUT.write('\n')

ID.close()
OUT.close()


