#!/usr/bin/python3

#-------------------------------------------------------
#-----------------merge the overlappinng peptide--------
#--------------------------------------------------------

def peptide_merge_lr(pep_list):
  #print('merge overlapping peptides to the longest ones')
  #print(pep_list)
  line_mgg=[]
  while(pep_list!=[]):
    pep_list2=[]
    pep_list3=[]
    corf=pep_list[0][-1]
    for jj in range(0,len(pep_list)):
      if (pep_list[jj][-1]==corf):
        pep_list2.append(pep_list[jj])
      else:
        pep_list3.append(pep_list[jj])
    pep_list=pep_list3
    line_mg=['','','','','','','']
    for jj in range(3,len(pep_list2[0])):
      line_mg[jj]=pep_list2[0][jj]

    ll=pep_list2[0][-1].find(pep_list2[0][1])
    rr=pep_list2[0][-1].find(pep_list2[0][1])+len(pep_list2[0][1])
    
    for line in pep_list2:
      for jj in range(0,3):
        line_mg[jj]=line_mg[jj]+','+line[jj]
    
      if(line[1] in line[-1][ll:rr]):
         continue
      elif(line[-1][ll:rr] in line[1]):
        ll=line[-1].find(line[1])
        rr=line[-1].find(line[1])+len(line[1])
       
      else:
        if(ll>line[-1].find(line[1])):
          ll=line[-1].find(line[1])
        if(rr<line[-1].find(line[1])+len(line[1])):
          rr=line[-1].find(line[1])+len(line[1])
        
          
    line_mg.insert(2,line_mg[-1][ll:rr])
    line_mgg.append(line_mg)
  return line_mgg

