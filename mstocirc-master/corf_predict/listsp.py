#!/usr/bin/python


def listsplit(sjj_list,step):
  sjj_list_list=[[]]
  i=0
  cun=0
  for sjj_line in sjj_list:
    if(cun>=step):
      i=i+1
      sjj_list_list.append([])
      cun=0
    sjj_list_list[i].append(sjj_line)
    cun=cun+1
  return sjj_list_list
