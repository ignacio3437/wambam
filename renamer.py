#!/usr/bin/python
import os,sys,re
filein=os.path.abspath("/Users/josec/Desktop/rbset2and3/rb2.tre")
filename=os.path.basename(filein)
pwd=os.path.dirname(filein)
rfilename="r_"+filename

codex="""W113582,CHANGE122284
W122213,CHANGE88644
W113636,CHANGE103233
W125322,CHANGE107942
W88713,CHANGE107977
W105891,CHANGE103907
W105887,CHANGE122241
W88704,CHANGE109117
W105892,CHANGE109131
W103202,CHANGE122277
W125321,CHANGE126259
W103263,CHANGE126280
W91917,CHANGE98774
W107988,CHANGE99667
W103216,CHANGE105213
W103271,CHANGE92455
W88646,CHANGE113636
W88700,CHANGE88704
W102824,CHANGE88708
W102827,CHANGE105892
W103200,CHANGE103210
W105890,CHANGE125321
W92127,CHANGE103202
W122287,CHANGE103263
W135688,CHANGE122218
W135656,CHANGE125322
W135677,CHANGE88698
W98771,CHANGE91917
W101156,CHANGE88713
W102826,CHANGE107988
W88653,CHANGE122213
W122218,CHANGE105891
W88708,CHANGE103216
W88698,CHANGE105887
W103210,CHANGE103271
W135658,CHANGE135658
W135659,CHANGE135659
W135671,CHANGE135671
W135672,CHANGE135672
W135673,CHANGE135673
W135721,CHANGE135721
W135687,CHANGE135687
W119987,CHANGE119987
W135720,CHANGE135720
W135678,CHANGE135678
W135679,CHANGE135679
W135693,CHANGE135693
W135681,CHANGE135681
W135682,CHANGE135682
W88641,CHANGE88641
W103942,CHANGE103942
W107184,CHANGE107184
W110368,CHANGE110368
W135514,CHANGE135514
W135519,CHANGE135519
W135592,CHANGE135592
W88705,CHANGE88705
W135566,CHANGE135566
W135567,CHANGE135567
W110363b,CHANGE110363b
W103963,CHANGE103963
W88710,CHANGE88710
W107968,CHANGE107968
W103899,CHANGE103899
W103948,CHANGE103948
W103908,CHANGE103908
W113566,CHANGE113566
W113576,CHANGE113576
W107182,CHANGE107182
W98875,CHANGE98875
W109135,CHANGE109135
W113654,CHANGE113654
W119973,CHANGE119973
W88714,CHANGE88714
W88706,CHANGE88706
W135676,CHANGE135676
W135657,CHANGE135657
W98765,CHANGE98765
W101166,CHANGE101166
W101165,CHANGE101165
W110363,CHANGE110363
W101157,CHANGE101157
W101164,CHANGE101164
W101160,CHANGE101160
W102154,CHANGE102154
W102155,CHANGE102155
W123301,CHANGE123301
W122895,CHANGE122895
W102832,CHANGE102832
W110359,CHANGE110359
W113574,CHANGE113574
W110356,CHANGE110356
W107982,CHANGE107982
W130374,CHANGE130374
W110366,CHANGE110366"""
clist=codex.split('\n')
codexdict={}
for c in clist:
    pair=c.split(',')
    codexdict[pair[0]]=pair[1]

with open(filein,'r') as filein:
    with open(os.path.join(pwd,rfilename),'w') as fileout:
        instring=filein.read()
        # fileout.write(instring)
        for key in codexdict:
            instring=instring.replace(key,codexdict[key])
        # fileout.write("\n\n\n")
        instring=instring.replace("CHANGE","W")
        fileout.write(instring)
print "done"