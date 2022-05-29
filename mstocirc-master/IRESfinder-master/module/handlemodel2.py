#1/usr/bin/python3

#adapted by 'qumu00'

import sys

finput = sys.argv[1]
foutput = sys.argv[2]
foutlab = sys.argv[3]
ws = int(sys.argv[4])
step = int(sys.argv[5])

#print('111'+finput+'\n')
ID = open(finput,'r+')
OUT = open(foutput,'w+')
LAB = open(foutlab,'w+')
idd=''
strr=''
for ID_line in ID.readlines():
  if(ID_line[0]=='>'):
    if(strr!=''):
      lenn =len(strr)
      if(lenn <= ws):
        LAB.write(idd+'\n')
        OUT.write(strr+'\n')  
      else:
        st = 0
        m = 0
        while(st + ws <lenn):
          wd = strr[st:st+ws]
          LAB.write(idd+'_')
          LAB.write('%d_'%(st+1))
          LAB.write('%d\n'%(st+ws))
          OUT.write('%s\n'%wd)
          st =st + step
          m=m+1
        st = lenn - ws
        wd = strr[st:st+ws]
        LAB.write('%s_'%idd)
        LAB.write('%d_'%(st+1))
        LAB.write('%d\n'%lenn)
        OUT.write(wd+'\n')
    
    strr = ''
    idd = ID_line
    idd=idd[1:].strip()
    
  else:
    strr=strr+ID_line
    strr=strr.strip()
      
if(strr!=''):
  lenn=len(strr)
  if(lenn <= ws ):
    LAB.write(idd+'\n')
    OUT.write(strr+'\n')
  else:
    st = 0
    m = 0
    while(st+ws < lenn):
      wd=strr[st:st+ws]
      LAB.write(idd+'_')
      LAB.write(str(st+1)+'_')
      LAB.write(str(st+ws)+'\n')
      OUT.write(wd+'\n')
      st=st+step
      m=m+1
    st = lenn - ws
    wd = strr[st:st+ws]
    LAB.write(idd+'_')
    LAB.write(str(st+1)+'_')
    LAB.write(str(lenn)+'\n')
    OUT.write(wd+'\n')



ID.close()
OUT.close()
LAB.close()
