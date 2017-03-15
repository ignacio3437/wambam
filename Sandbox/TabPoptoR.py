#!/usr/bin/python


Pops="""StraitsMagellan
BurdwoodBanks
FalklandIsl
ShagRocks
SouthGeorgia
SouthSandwich
HerdmanBank
DiscoveryBank
SoutOrkney
ElephantIsl1
ElephantIsl2
BransStrait1
BransStrait2
ScottA
RossSea
PrydzBay
ShelfBreak
HeardIslCo
HeardIslAu"""
with open("/Users/josec/Desktop/Porania/PoraniaSubset2/PoraniaSubset.tsv",'rU') as tsvh:
	tsv=tsvh.readlines()
popdict={}
for line in tsv:
	linel=line.strip('\n').split('\t')
	try:
		popdict[linel[1]]
		popdict[linel[1]].append(linel[0])
	except:
		popdict[linel[1]]=[linel[0]]
for x in Pops.split('\n'):
	namestring='","'.join(popdict[x])
	print '%s <- c("%s")'%(x,namestring)
popstring=','.join(Pops.split('\n'))
print "GENOME.class <- set.populations(GENOME.class,list(%s))"%(popstring)
# pop.1 <- c("seq1","seq2")

"""
GENOME.class <- readData("/Users/josec/Desktop/Fasta")
StraitsMagellan <- c("S20019","S20020","S20021","S20024","S20071","S20072","S20073","S20076","S20077","S20078","S20141","S20142","S20143","S20144","S20145","S20146","S20147","S20148","S20149","S20150")
BurdwoodBanks <- c("S5576","S5577","S5578","S5579","S5580","S5581","S5582","S5583","S5584","S5585","S5586","S5636","S5637","S5638","S5639","S5640","S5641","S5642","S5643","S5644","S5645","S5646","S5647","S5648","S5649","S5650","S5651","S5652","S5653","S5654","S5655")
FalklandIsl <- c("S20839","S20840","S20871","S20872","S20873","S20874","S20925","S20926","S20927","S20929","S20930","S20931","S20932","S20933","S20934","S20939","S20940","S20978")
ShagRocks <- c("S20364","S20365","S20371","S20372","S20405","S20406","S20426","S20430","S20476")
SouthGeorgia <- c("S3284","S3285","S3286","S3287","S3288","S3289","S3290","S3291","S3912","S3913","S3914","S3915","S3916","S3917","S3918","S3919","S3920","S3921","S3922","S3923","S3924","S3925","S3926","S3928","S3929")
SouthSandwich <- c("S0168","S0169","S0170","S0171","S0476","S0477","S0478","S0519","S0520","S0335","S0336","S0597","S0798")
HerdmanBank <- c("S4401","S4402","S4403","S4404","S4421","S4425","S4426","S4428","S4485")
DiscoveryBank <- c("S1676","S1679","S1837","S1838","S1839","S1852","S1853","S1978","S1980","S4012","S4106","S4107","S4114","S4115","S4116","S4170","S4301","S4302","S4303")
SoutOrkney <- c("S1666","SO1004","A","B","C","D","E","F","G","H","I","J","K","L","M")
ElephantIsl1 <- c("S6871","S6872","S6873","S6874","S6875","S6876","S6877","S6878","S6879","S6880","S6881","S6882","S6883")
ElephantIsl2 <- c("S6601","S6602","S6603","S6604","S6605","S6606","S6608","S6609","S6610","S6611","S6612","S6613","S6614","S6615","S6616","S6617","S6618","S6619","S6620")
BransStrait1 <- c("S4129","S4130","S4131","S4132","S4133","S4134","S4702","S4703","S4704","S6686","S6687","S6688","S4705")
BransStrait2 <- c("S5956","S5958","S5959","S5960","S5961","S5962","S5964","S5965","S5966","S5967","S5968","S5969","S5970","S6293","S6294","S6295","S6296","S6298","S6299")
ScottA <- c("N0044","N0045","N0046","N0047","N0048","N0022","N0029","N0030","N0013","N0026","N0027","N0028","N0031","N0032","N0033","N0049","N0050","N0051","N0052","N0053","N0054","N0055","N0056","N0057","N0058","N0059","N0060","N0061","N0062","N0063","N0064","N0021")
RossSea <- c("N0034","N0035","N0036","N0037","N0038","N0039","N0040","N0041","N0042","N0043","N0023","N0024","N0025")
PrydzBay <- c("AA113","AA114","AA115","AA116","AA117","AA118","AA119","AA120","AA121","AA124","AA125","AA126","AA127","AA128","AA130","AA131","AA132","AA133","AA159","AA160")
ShelfBreak <- c("PNG540","PNG541","PNG542","PNG543","PNG544","PNG532","PNG533","PNG534","PNG536","PNG537","PNG538","PNG539","PNG535","PNG545","PNG546","PNG547","PNG548","PNG549","PNG550","PNG551")
HeardIslCo <- c("PNG559","PNG560","PNG561","PNG557","PNG558","AAJ097")
HeardIslAu <- c("AAJ080","AAJ083","AAJ084","AAJ085","AAJ086","AAJ087","AAJ088","AAJ089","AAJ090","AAJ091","AAJ095","AAJ096","PNG556")
GENOME.class <- set.populations(GENOME.class,list(StraitsMagellan,BurdwoodBanks,FalklandIsl,ShagRocks,SouthGeorgia,SouthSandwich,HerdmanBank,DiscoveryBank,SoutOrkney,ElephantIsl1,ElephantIsl2,BransStrait1,BransStrait2,ScottA,RossSea,PrydzBay,ShelfBreak,HeardIslCo,HeardIslAu))
GENOME.class <- F_ST.stats(GENOME.class,list(StraitsMagellan,BurdwoodBanks,FalklandIsl,ShagRocks,SouthGeorgia,SouthSandwich,HerdmanBank,DiscoveryBank,SoutOrkney,ElephantIsl1,ElephantIsl2,BransStrait1,BransStrait2,ScottA,RossSea,PrydzBay,ShelfBreak,HeardIslCo,HeardIslAu))
GENOME.class@nuc.F_ST.pairwise
"""