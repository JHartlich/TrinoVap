#!/usr/bin/env python
# -*- coding: utf-8 -*-

#---packages----------------------
import argparse, os, sys

#---functions---------------------
## creating a command line
def command():
  Parser = argparse.ArgumentParser(prog='TrinoVap.py', usage='./%(prog)s TransDecoder FASTA Trinotate annotation report',
    description='%(prog)s is a simple tool to check for SwissProt, Kegg or Pfam annotations for a transcript based on the Trinotate annotation report.\n%(prog)s uses any TransDecoder output in FASTA format to report annotions for complete, internal and partial sequences seperately.',
    formatter_class=argparse.RawTextHelpFormatter)
  Parser.add_argument('INPUT', nargs=2, help='name of TransDecoder FASTA (must be first argument), \nname of Trinotate annotation report TXT file (must be second argument), \nplease add directory if file is not in current working directory')
  Args = Parser.parse_args()
  return (Args)

## checking for SwissProt annotation
def anno_sprot(sphit):
  if sphit !='.':
    if sphit.find('`') != -1:
      sp = sphit.split('`')[0]
      e_val = sp.split('^')[4]
      if e_val.find('e-') != -1:
        value = e_val.split('-')[1]
        if int(value) >= 5:
          spname = sp.split('^')[5]
          spname = spname.split('=')[1]
          return ('%s;%s' % (spname.strip(';'),e_val.strip('E:')))
        else:
          return ('.')
      else:
        return ('.')
    else:
      e_val = sphit.split('^')[4]
      if e_val.find('e-') != -1:
        value = e_val.split('-')[1]
        if int(value) >= 5:
          spname = sphit.split('^')[5]
          spname = spname.split('=')[1]
          return ('%s;%s' % (spname.strip(';'),e_val.strip('E:')))
        else:
          return ('.')
      else:
        return ('.')
  else:
    return (sphit)

## checking for Kegg annotation
def anno_kegg(khit):
  if khit.find('KO') != -1:
    return (khit.split('KO:')[1])
  else:
    return ('.')

## checking for Pfam annotation
def anno_pfam(pfhit):
  if pfhit != '.':
    if pfhit.find('`') != -1:
      pf = pfhit.split('`')[0]
      e_val = pf.split('^')[4]
      if e_val.find('e-') != -1:
        value = e_val.split('-')[1]
        if int(value) >= 5:
          pfname = pf.split('^')[2]
          return ('%s;%s' % (pfname,e_val.strip('E:')))
        else:
          return ('.')
      else:
        return ('.')
    else:
      e_val = pfhit.split('^')[4]
      if e_val.find('e-') != -1:
        value = e_val.split('-')[1]
        if int(value) >= 5:
          pfname = pfhit.split('^')[2]
          return ('%s;%s' % (pfname,e_val.strip('E:')))
        else:
          return ('.')
      else:
        return ('.')
  else:
    return (pfhit)


#---main program------------------
## checking command line
ARGs = command()

## generating outputfiles
NewCompl = 'Trinotate_complete_peps.tsv'
NewC = open(NewCompl,'w')
NewC.write('Transcript\tSwissProt\tKegg\tPfam')

NewInter = 'Trinotate_internal_peps.tsv'
NewI = open(NewInter,'w')
NewI.write('Transcript\tSwissProt\tKegg\tPfam')

New5prim = 'Trinotate_fiveprimepart_peps.tsv'
New5 = open(New5prim,'w')
New5.write('Transcript\tSwissProt\tKegg\tPfam')

New3prim = 'Trinotate_threeprimepart_peps.tsv'
New3 = open(New3prim,'w')
New3.write('Transcript\tSwissProt\tKegg\tPfam')

NewStats = 'TrinoVap_stats.txt'
NewSt = open(NewStats, 'w')

## accessing FASTA file
Ffile = open(ARGs.INPUT[0],'r')
FAA = Ffile.read()
Ffile.close()
FAA = FAA.replace('\r','')
FAA = FAA.split('>')[1:] # seperating sequences -> list

Ctype = []
Itype = []
Ftype = []
Ttype = []
### extracting sequences
for entry in FAA:
  Label = entry.split('\n')[0]
  Gene = Label.split(' ')[0]
  if Label.find('type:comp')!=-1:
    Ctype[len(Ctype):] = [Gene]
  elif Label.find('type:int')!=-1:
    Itype[len(Itype):] = [Gene]
  elif Label.find('type:5')!=-1:
    Ftype[len(Ftype):] = [Gene]
  elif Label.find('type:3')!=-1:
    Ttype[len(Ttype):] = [Gene]

FAA = ''

## accessing Trinotate file
File = open(ARGs.INPUT[1],'r')
TRF = File.read()
File.close()
TRF = TRF.replace('\r','')
TRF = TRF.strip('\n')
#print(TRF.split('\n')[0])
TRF = TRF.split('\n')[1:]

## initiating counter for statistics
sp_count_c = 0
sp_count_i = 0
sp_count_5 = 0 
sp_count_3 = 0
pf_count_c = 0 
pf_count_i = 0
pf_count_5 = 0
pf_count_3 = 0
kg_count_c = 0
kg_count_i = 0
kg_count_5 = 0
kg_count_3 = 0

## extrakt information!
for line in TRF:
  gene = line.split('\t')[4]
  # sprot = line.split('\t')[6] # blastp hits for SwissProt
  sprot = line.split('\t')[2] # blastx hits for SwissProt
  pfam  = line.split('\t')[7]
  kegg  = line.split('\t')[11]

  if gene in Ctype:
    Sprot = anno_sprot(sprot)
    if Sprot != '.':
      sp_count_c += 1
      NewC.write('\n%s\t%s\t.\t.' % (gene,Sprot))
    else:
      Kegg = anno_kegg(kegg)
      if Kegg != '.':
        kg_count_c += 1
        NewC.write('\n%s\t.\t%s\t.' % (gene,Kegg))
      else:
        Pfam = anno_pfam(pfam)
        if Pfam != '.':
          pf_count_c += 1
          NewC.write('\n%s\t.\t.\t%s' % (gene,Pfam))

  elif gene in Itype:
    Sprot = anno_sprot(sprot)
    if Sprot != '.':
      sp_count_i += 1
      NewI.write('\n%s\t%s\t.\t.' % (gene,Sprot))
    else:
      Kegg = anno_kegg(kegg)
      if Kegg != '.':
        kg_count_i += 1
        NewI.write('\n%s\t.\t%s\t.' % (gene,Kegg))
      else:
        Pfam = anno_pfam(pfam)
        if Pfam != '.':
          pf_count_i += 1
          NewI.write('\n%s\t.\t.\t%s' % (gene,Pfam))

  elif gene in Ftype:
    Sprot = anno_sprot(sprot)
    if Sprot != '.':
      sp_count_5 += 1
      New5.write('\n%s\t%s\t.\t.' % (gene,Sprot))
    else:
      Kegg = anno_kegg(kegg)
      if Kegg != '.':
        kg_count_5 += 1
        New5.write('\n%s\t.\t%s\t.' % (gene,Kegg))
      else:
        Pfam = anno_pfam(pfam)
        if Pfam != '.':
          pf_count_5 += 1
          New5.write('\n%s\t.\t.\t%s' % (gene,Pfam))

  elif gene in Ttype:
    Sprot = anno_sprot(sprot)
    if Sprot != '.':
      sp_count_3 += 1
      New3.write('\n%s\t%s\t.\t.' % (gene,Sprot))
    else:
      Kegg = anno_kegg(kegg)
      if Kegg != '.':
        kg_count_3 += 1
        New3.write('\n%s\t.\t%s\t.' % (gene,Kegg))
      else:
        Pfam = anno_pfam(pfam)
        if Pfam != '.':
          pf_count_3 += 1
          New3.write('\n%s\t.\t.\t%s' % (gene,Pfam))

NewSt.write('\nCOMPLETE    :\tSprot: %d\tPfam: %d\tKegg: %d' %(sp_count_c,pf_count_c,kg_count_c))
NewSt.write('\nINTERNAL    :\tSprot: %d\tPfam: %d\tKegg: %d' %(sp_count_i,pf_count_i,kg_count_i))
NewSt.write('\n5 PRIME PART:\tSprot: %d\tPfam: %d\tKegg: %d' %(sp_count_5,pf_count_5,kg_count_5))
NewSt.write('\n3 PRIME PART:\tSprot: %d\tPfam: %d\tKegg: %d' %(sp_count_3,pf_count_3,kg_count_3))

NewC.close()
NewI.close()
New5.close()
New3.close()
NewSt.close()