#!/usr/bin/python3

def proano(file_annote,curpath):
  print('extract the exon positon from the annotion.gtf...')
  #fi=open(args.annotion,'r+')
  fi=open(file_annote,'r+')
  fo=open(curpath+'/exon_infor.txt','w+')
  fp=open(curpath+'/exon_infor2.txt','w+')
  lock_gene=0
  llock=0
  exon_cstr=[]
  exon_cend=[]
  exon_str2=[]
  exon_end2=[]
  chrm=[]
  chrm2=[]
  std=''
  stdl=[]
  stdl2=[]
  cuunt=0
  #for sp in fi.readlines():
  while(1):
    sp=fi.readline()
    if(sp==''):
      break

    sp=sp.strip()
    spp=sp.split('\t')
    #print(spp)
    if(sp[:2]=='#!' or sp[0]=='#'): 
      continue 
    #print '03'
    fp.write('chr'+spp[0]+'\t'+spp[2]+'\t'+spp[3]+'\t'+spp[4]+'\t%s\n'%spp[6])
    if(spp[2]=='gene'):
      lock_gene=lock_gene+1
    if(spp[2]=='transcript' or spp[2]=='mRNA' or spp[2]=='lnc_RNA'):
      #print('transcript')
      llock=1
    if(lock_gene%2==0):
      llock=1
      lock_gene=lock_gene-1
      
    if(spp[2]=='exon' and llock==0):
      #print(spp[3],spp[4])
      chrm.append('chr'+spp[0])
      chrm2.append('chr'+spp[0])
      exon_str2.append(spp[3])
      exon_end2.append(spp[4])
      #std=spp[6]
      stdl.append(spp[6])
      stdl2.append(spp[6])
      
      
    if(llock==1):
      num=len(exon_str2)
      llock=0
      
      if(num==0):
        continue
      for i in range(0,num):
          
        if(stdl[i]=='+'):
          exon_cstr.append(exon_str2[i])
          exon_cend.append(exon_end2[i])
          #print '01'
          fo.write(chrm[i]+'\texon\t'+exon_str2[i]+'\t'+exon_end2[i]+'\t'+stdl[i]+'\n')
        else:
          exon_cstr.append(exon_str2[num-1-i])
          exon_cend.append(exon_end2[num-1-i])
          fo.write(chrm[num-1-i]+'\texon\t'+exon_str2[num-1-i]+'\t'+exon_end2[num-1-i]+'\t'+stdl[num-1-i]+'\n')
      
      exon_str2=[]
      exon_end2=[]
      chrm=[]
      stdl=[]
      #chrm.append('~')
      chrm2.append('~')
      #stdl.append('~')
      stdl2.append('~')
      exon_cstr.append('~')
      exon_cend.append('~')
      fo.write('~\texon\t~\t~\t~\n')
    
  num=len(exon_str2)
 
  for i in range(0,num):
    if(stdl[i]=='+'):
      exon_cstr.append(exon_str2[i])
      exon_cend.append(exon_end2[i])
          #print '01'
      fo.write(chrm[i]+'\texon\t'+exon_str2[i]+'\t'+exon_end2[i]+'\t'+stdl[i]+'\n')
    else:
      exon_cstr.append(exon_str2[num-1-i])
      exon_cend.append(exon_end2[num-1-i])
      fo.write(chrm[num-1-i]+'\texon\t'+exon_str2[num-1-i]+'\t'+exon_end2[num-1-i]+'\t'+stdl[num-1-i]+'\n')
  #print(len(chrm2),len(exon_cstr),len(exon_cend),len(stdl2))
  #for i in range(0,len(exon_cstr)):
  #  fo.write(chrm[i]+'\t'+'exon'+'\t'+exon_cstr[i]+'\t'+exon_cend[i]+'\t'+stdl[i]+'\n')
  fi.close()
  fo.close()
  fp.close()
  return chrm2,exon_cstr,exon_cend,stdl2
