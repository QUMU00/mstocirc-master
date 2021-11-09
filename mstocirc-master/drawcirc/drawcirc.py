#!/usr/bin/python3
__Author__='QUMU00'

import numpy
import matplotlib.pyplot as pltt
import random 
import sys

# b---blue  c---cyan  g---green  k---black
# m---magenta r---red  w---white yellow
# tamato  sanddlebrown  darkorange  gold
# palegreen lurquoise deepskyblue  darkviolet
#
# https://finthon.com/matplotlib-color-list/

listcolor=['red','cyan','green','black','yellow','gold','darkviolet']

a,b=0.,0.

def line(ang,rr):
  rra=numpy.arange(rr-0.25,rr+0.25,0.01)
  lx=[]
  ly=[]
  anga=numpy.arange(ang-0.03,ang+0.03,0.01)
  #print(anga)
  for n in range(0,len(anga)-1):
    for m in range(0,len(rra)-1):
      lx.append(a+rra[m]*numpy.sin(anga[n]))
      ly.append(b+rra[m]*numpy.cos(anga[n]))
  
  return lx,ly


def drawcircle(strr,px):
  
  #print("plese input the parameter of circRNA\nlike this:")
  #print("[circ_001385,678,56,921,125,225,0-341-451-537-678,678nt expected prorein 195aa 17KD]")
  #nam='circ_001385'
  #info='lengnt\nexpected prorein 195aa\n17KD'
  #leng=678
  #orf_st=56
  #orf_ed=921
  #ires_s=125
  #ires_e=225
  
  #strr=input(":")
  #print(strr)
  #pltt.plot(211)
  strr=strr.split(',')[:8]
  #print(strr)
  nam=strr[0]
  info=strr[7]
  leng=int(strr[1])
  orf_st=int(strr[2])
  orf_ed=int(strr[3])
  ires_s=int(strr[4])
  ires_e=int(strr[5])
  
  r=numpy.arange(1.9,2.1,0.01)
  wr=2.5
  
  fig=pltt.figure()
  axes=fig.add_subplot(111)
  axes.set_xlim(-3,3)
  axes.set_ylim(-3,3)

  #draw exon position
  #exon='0,341,451,537,678'
  exon=strr[6]
  exonl=exon.split('-')
  for m in range(0,len(exonl)-1):
    theta=numpy.arange(int(exonl[m])/leng*2*numpy.pi,int(exonl[m+1])/leng*2*numpy.pi,0.01)
    for i in range(0,len(r)):
      x=a+r[i]*numpy.sin(theta)
      y=b+r[i]*numpy.cos(theta)
      axes.plot(x,y,listcolor[m%7])
  
  x,y=line(0,2.0)
  axes.plot(x,y,'yellow')
  pltt.text(x[0]+0.15,y[0],'junction site')
  
  #draw  open reading frame,orf
  wb=b-0.1
  wr=wr-0.1
  orf_s=orf_st/leng
  orf_e=orf_ed/leng
  x,y=line(orf_s*2*numpy.pi,2.0)
  axes.plot(x,y,'green')
  pltt.text(x[len(x)-1]+0.15,y[len(y)-1],'ATG,initation translation site')
  x,y=line(orf_e*2*numpy.pi,2.0)
  axes.plot(x,y,'red')
  pltt.text(x[len(x)-1]+0.15,y[len(y)-1],'TGA,translation termation site')
  
  
  orf_s2=orf_s
  j=1
  theta2=numpy.arange(0,2,0.1)
  xx_x,xx_y=0,0
  while(orf_s2<=orf_e):
    if(orf_s2<=0.5*j):
      if(orf_s2+0.5>=orf_e):
        theta2=numpy.arange(orf_s2*2*numpy.pi,orf_e*2*numpy.pi,0.01)
      else:
        theta2=numpy.arange(orf_s2*2*numpy.pi,j*numpy.pi,0.01)
      wr=wr+0.1
      if(j%2==1):
        wb=wb-0.1
      else:
        wb=wb+0.1
      x=a+wr*numpy.sin(theta2)
      y=wb+wr*numpy.cos(theta2)
      axes.plot(x,y,'g')
      if(xx_y==0):
         xx_x,xx_y=x[0],y[0]
      orf_s2=0.5*j
    j=j+1
  #jj=random.randint(0,len(x)-1)
  pltt.text(xx_x+0.05,xx_y,'open reading frame')

  # draw Ires
  ir=2.25
  theta3=numpy.arange(ires_s*2/leng*numpy.pi,ires_e*2/leng*numpy.pi,0.01)
  x=a+ir*numpy.sin(theta3)
  y=b+ir*numpy.cos(theta3)
  axes.plot(x,y,'m')
  pltt.text(x[0]+0.05,y[0],'+%d'%ires_s)
  jj=random.randint(0,len(theta3)-1)
  pltt.text(x[jj]+0.05,y[jj],'IRES')
  pltt.text(x[len(theta3)-1]+0.05,y[len(theta3)-1],'+%d'%ires_e)

  x,y=line(ires_s/leng*2*numpy.pi,2.0)
  axes.plot(x,y,'darkviolet')
  x,y=line(ires_e/leng*2*numpy.pi,2.0)
  axes.plot(x,y,'darkviolet')
  
  #draw circRNA name
  pltt.text(-0.5,-3,nam)
  pltt.text(-0.8,0.3,info)
  
  axes.axis('equal')
  axes.axis('scaled')
  pltt.savefig(px+'/'+nam+'_circc.png')
  #pltt.show()


#---------------------------------------------------------------
#-----------------------------------------------------------------
#---------------------------------------------------------------

import math



def drawcirclew(strr,px):     
  strr=strr.split(',')
  #dd=input('Under performance testing,any button to continue:')
  #print('-------------------------------------')
  #ag=['A','G','C','T']
  circ=strr[8]
  ll=len(circ)
  cORF=strr[9].split('-')[1]
  pep=strr[9].split('-')[0]
  pep_ind=cORF.find(pep)*3+int(strr[2])
  #print(pep_ind[1],int(strr[1]))
  #nm=list(map(int,strr[6].split('-')))
  #nm=[0,int(strr[3])-int(strr[1]),int(strr[4]),int(strr[5]),int(strr[2]),int(strr[1])]
  nm=[0,int(strr[3])-int(strr[1]),int(strr[2]),int(strr[1])]
  #print(nm)
  if(nm[1]>=nm[2]):
    nm[2]=nm[1]
    cORF=cORF[:ll//3-3]+'...'
  #aa=['A','R','D','C','Q','E','H','I','G','N','L','K','M','F','P','S','T','W','Y','V']
  #for i in range(0,663):
  #  circ=circ+ag[random.randint(0,3)]
  #nm=[1,50,127,200,350,663]
  #    circ_start corf_end ires_start ires_end corf_start circ_end
  #for i in  range(0,(663+50-350)//3-1):
  #  cORF=cORF+aa[random.randint(0,19)]
  #print(cORF)
  # test staiistcics ,to randomly generate the dna sequence anf peptide sequence of proper length 
  #print(nm) 
  rp=[]
  rp2=[]
  print('orgn',ll)
  if(ll*0.8>280):
    print('1')
    for i in range(0,len(nm)-1):
      if(nm[i+1]-nm[i]>ll*0.2):  
        li=(ll-280)*(nm[i+1]-nm[i]+1)//ll+3
        #print(nm[i],nm[i+1]-li-2)
        rr=random.randint(nm[i],nm[i+1]-li-2)
        #rp.append(circ[rr:rr+li])
        rp.append(rr)
        rp2.append(li)
      else:
        rp.append(0)
        rp2.append(0)
    lrp=0  
    for i in range(0,len(rp)):
      if(rp[i]!=0):
        rp[i]=rp[i]-lrp
        #circ=circ.replace('...',rp[i])
        print(rp[i],circ)
        circ=circ[:rp[i]]+'...'+circ[rp[i]+rp2[i]:]
        lrp=rp2[i]+lrp-3
        nm[i+1]=nm[i+1]-lrp
        
    
    if(rp[0]!=0):
      cORF=cORF[:(int(strr[1])-int(strr[2])+rp[0])//3]+'...'+cORF[(int(strr[1])-int(strr[2])+rp[0]+rp2[0])//3:]
    if(rp[2]!=0):
      cORF=cORF[:(rp[2]-nm[2])//3]+'...'+cORF[(rp[2]+rp2[2]-nm[2])//3:]
    #rmm=random.randint(1,len(cORF)-len(rp[0]+rp[-1])//3-2)
    #cORF=cORF.replace(cORF[rmm:rmm+len(rp[0]+rp[-1])//3],'...')
  #cORF='012345678901234567890123456789012345'
  # to randomly delete partial sequences of aimed dna and peptide  
  ll=len(circ)
  print('del,',ll)
  fig = pltt.figure(figsize=(10, 10),dpi=80)
  ax = fig.add_subplot(111)
  ax.axis('equal')
  #pltt.ylim(-4,4)
  #pltt.xlim(-4,4)
  ax.set(xlim=[-6.0, 6.0], ylim=[-6.0, 6.0], title='circRNA base circle', ylabel='Y-Axis', xlabel='X-Axis',)
  
  #pltt.subplot(2,1,1)
  theta=numpy.linspace(0*numpy.pi,2.0*numpy.pi,ll+1)
  r=4.5
  a,b=0.,0.
  x=a+r*numpy.sin(theta)
  y=b+r*numpy.cos(theta)
  cl='black'
  #print(nm)
  ax.text(x[0],y[0]+0.2,'|',color='red')
  for i in range(0,ll):
    #yellow,red,black,green,blue,violet
    if(i<=nm[1] or i>=nm[-2]):
      cl='green'
    #if(nm[2]<=i and i<=nm[3]):
    #  cl='violet'
    ax.text(x[i],y[i],circ[i],fontsize=6,rotation=-360.0*i/ll,color=cl)
    cl='black'
  # draw the base sequence circle,green color for cORF and violet for IRES eement   


  theta=numpy.linspace(nm[-2]/nm[-1]*2*numpy.pi,(nm[-1]+nm[1])/nm[-1]*2*numpy.pi,len(cORF))
  ang=numpy.linspace(-360*nm[-2]/nm[-1],-360*(nm[-1]+nm[1])/nm[-1],len(cORF)) 
  xa=a+(r-1.5)*numpy.sin(theta)
  ya=b+(r-1.5)*numpy.cos(theta)
  for i in range(0,len(cORF)):
    ax.text(xa[i],ya[i],cORF[i],fontsize=10,rotation=ang[i])
  #draw the peptide circe,inside of the base circe
  
  #if(cORF.find(pep)>=0):
  #  pep_ind=cORF.find(pep)*3+nm[-2]
  #else:
  #  pep_ind=nm[-1]*pep_ind/int(strr[1])
  pep_ind=nm[-1]*pep_ind/int(strr[1])
  theta=numpy.linspace(pep_ind/nm[-1]*2*numpy.pi,(pep_ind+len(pep)*3)/nm[-1]*2*numpy.pi,len(pep))
  ang=numpy.linspace(-360*(pep_ind)/nm[-1],-360*(pep_ind+len(pep)*3)/nm[-1],len(pep)) 
  xa=a+(r-1.2)*numpy.sin(theta)
  ya=b+(r-1.2)*numpy.cos(theta)
  for i in range(0,len(pep)):
    ax.text(xa[i],ya[i],pep[i],fontsize=10,rotation=ang[i],color='blue')
  
  pltt.savefig(px+'/'+strr[0]+'_circw .png')
  #pltt.show()




#----------------------------------------------
#---------------main function------------------
#-----------------------------------------------

fi=open(sys.argv[1],'r+')
#cc=input('whether or not show the drawing picture:')
cc='no'
count=1
for si in fi.readlines():
  if(si[0]=='#'):
    continue
  if(int(si.split(',')[1])>1500 or int(si.split(',')[5])<50):
    continue
  drawcircle(si.strip(),sys.argv[2])
  drawcirclew(si.strip(),sys.argv[2])
  if(cc[0]=='y' or cc[0]=='Y'):
    pltt.show()
  if(count>=10):
    print('there will consume too many memory while more than 10 figurs is opened .\n\
you can set the input file to draw the specific circRNA.\nthankyou.')
    break
  count=count+1
fi.close()

