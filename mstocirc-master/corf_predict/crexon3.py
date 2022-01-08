#!/usr/bin/python3

def circRNA_extractp(exoo,fo_temp_list,fq_temp_list,list_list_sj,chrm_exon_dict):
  for line in list_list_sj:
    circRNA_extract(exoo,fo_temp_list,fq_temp_list,line,chrm_exon_dict)

#----------------------------------------
#-----------------------------------------------------------------------

def circRNA_extract(exoo,fo_temp_list,fq_temp_list,sjj_string,chrm_exon_dict):
  strddd=''
  sjt=['','','','','','','']
  sjtt=sjj_string.split('\t')
  print(sjt[0])
  str_ind=0
  exos_str=''
  exos_end=''
  #fo_temp_line=''
  #fq_temp_line=''
  #exoo_line=''

  #sjt: circID, null,chrm.start,end,gene_symbl,null
  #input: chrom start end strand circRNA_ID transcript gene_symbol
  sjt[0]=sjtt[4]
  sjt[2]=sjtt[0]
  sjt[3]=sjtt[1]
  sjt[4]=sjtt[2]
  sjt[5]=sjtt[6]
    
    
  for i in range(0,len(chrm_exon_dict['chr%s_exon_str'%sjt[2]])):
    if(chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]==sjt[3]):
      str_ind=i
      strddd=chrm_exon_dict['chr%s_strd'%sjt[2]][i]
      #fq.write(sjt[0]+'\t'+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]+'\t'+strddd+'\n')
      fq_temp_list.append(sjt[0]+'\t'+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]+'\t'+strddd+'\n')
      print('dfg-1')
      break
  if(str_ind==0):
      #continue
      return ['','','']
  #for i in range(str_ind,len(chrm_exon_dict['chr%s_exon_end'%sjt[2]])):
  for i in range(str_ind,str_ind+50):
    #print(sjt[4])
    if(i>=len(chrm_exon_dict['chr%s_exon_end'%sjt[2]]) or chrm_exon_dict['chr%s_exon_end'%sjt[2]][i]=='~'):
      return ['','','']
    #if():
    #  break
    if(chrm_exon_dict['chr%s_exon_end'%sjt[2]][i]==sjt[4]):
      print('dfg-2')
      exos_str=exos_str+','+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]
      exos_end=exos_end+','+sjt[4]
      exoo.append(sjt[0]+' '+'chr'+sjt[2]+' '+sjt[3]+' '+sjt[4]+' %s '%sjt[5]+strddd+' '+exos_str[1:]+' '+exos_end[1:])
      #exoo_line=sjt[0]+' '+'chr'+sjt[2]+' '+sjt[3]+' '+sjt[4]+' %s '%sjt[5]+strddd+' '+exos_str[1:]+' '+exos_end[1:]
      #fo.write(sjt[0]+' '+'chr'+sjt[2]+' '+sjt[3]+' '+sjt[4]+' %s '%sjt[5]+strddd+' '+exos_str[1:]+' '+exos_end[1:]+'\n')
      fo_temp_list.append(sjt[0]+' '+'chr'+sjt[2]+' '+sjt[3]+' '+sjt[4]+' %s '%sjt[5]+strddd+' '+exos_str[1:]+' '+exos_end[1:]+'\n')
      break
      
    if(chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]==sjt[3]):
      exos_str=''
      exos_end=''
    exos_str=exos_str+','+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]
    exos_end=exos_end+','+chrm_exon_dict['chr%s_exon_end'%sjt[2]][i]
 
  #return [fo_temp_line,fq_temp_line,exoo_line]
  #return fo_temp_line,fq_temp_line,exoo_line
   
