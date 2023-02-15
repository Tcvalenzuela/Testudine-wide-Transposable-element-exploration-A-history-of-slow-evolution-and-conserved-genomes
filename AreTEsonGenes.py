import argparse
parser = argparse.ArgumentParser(description= "Identify if TE are in genetic feature and frag the kind of interaction")
parser.add_argument("--TEbed", "-bedfileforTE",  help="bed file with the information of TE")
parser.add_argument("--gff", "-gfffile",  help="Gff file with the genetic features")
#parser.add_argument("--O", "-output_file",  help="File to write the output to")
arg = parser.parse_args()
TEbed = open(arg.TEbed)
gff = open(arg.gff)
TEStarts={}
TEEnds={}
TEdNdS={}
TEsize={}
#print ("reading TEs on "+ str(arg.TEbed))
for TE in TEbed:
    if TE.startswith("Scaffold"):
        continue
    Splitted=TE.rstrip("\n").split("\t")
    Scaff=Splitted[0].split(".")[0]
    Start=int(Splitted[1])
    End=int(Splitted[2])
    Kind=str(Splitted[3])
    Kimura=float(Splitted[4])
    dNdS=float(Splitted[5])
    Size=End-Start
    TEStarts[TE]=(Start,Scaff)
    TEEnds[TE]=(End,Scaff)
    TEdNdS[TE]=(dNdS,Kimura)
    TEsize[TE]=(Start,End)
Gffexon={}
GffCDS={}
Gffgene={}
GffmRNA={}
Gfftranscript={}
GffcDNA_match={}
#print ("Done reading TEs, now readinf Gff "+ str(arg.gff) )
for line in gff:
    if line.startswith("#"):
        continue
    Split=line.split("\t")
    gffScaff=Split[0].split(".")[0]
    gffclass=Split[2]
    gffStart=int(Split[3])
    gffEnd=int(Split[4])
    if gffclass == "region":
        continue
    elif gffclass == "exon":
        Gffexon[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "CDS":
        GffCDS[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "mRNA":
        GffmRNA[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "gene":
        Gffgene[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "transcript":
        Gfftranscript[line]=(gffStart,gffEnd,gffScaff)
    elif gffclass == "cDNA_match":
        GffcDNA_match[line]=(gffStart,gffEnd,gffScaff)
print ("TE_scaff" +"\t"+ "TE_Start"+"\t"+ "TE_Ends" +"\t"+ "TE_familie" +"\t"+ "Kimura" +"\t"+ "DnDs" +"\t"+ "Mutation_rate" +"\t"+ "Kind of interaction"+"\t" +"gene_Scaff" + "\t" +"Gff_class" +  "\t" +"Gene_Starts" +"\t" +"Gene_Ends"+"\t" +"Score"+"\t" +"Strand"+"\t" +"phase"+"\t" +"Gene_ID_info")
TEsinGENES={}
for ele in Gffgene:
    for rec in TEStarts:
        if Gffgene[ele][2]!= TEStarts[rec][1]:
            continue
        elif Gffgene[ele][0] > TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec][0] and  Gffgene[ele][1] < TEEnds[rec][0] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "gene" + "_inside_TE"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("gene_inside_TE")
        elif Gffgene[ele][0] < TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec][0] and  Gffgene[ele][1] > TEEnds[rec][0] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "\t"+"TE_inside_"+ "gene" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_inside_gene")
        elif Gffgene[ele][0] <= TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec][0] and  Gffgene[ele][1] <= TEEnds[rec][0] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2Downstream_"+ "gene" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2Downstream_gene")
        elif Gffgene[ele][0] >= TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec][0] and  Gffgene[ele][1] >= TEEnds[rec][0] and TEStarts[rec][0] < Gffgene[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2upstream_"+"gene" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2upstream_gene")
        elif Gffgene[ele][0] == TEStarts[rec][0] and Gffgene[ele][0] < TEEnds[rec][0] and  Gffgene[ele][1] == TEEnds[rec][0] and TEStarts[rec][0] < Gffgene[ele][1] :
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_exactly_"+ "gene"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_exactly_gene")

for ele in GffmRNA:
    for rec in TEStarts:
        if GffmRNA[ele][2]!= TEStarts[rec][1]:
            continue
        elif GffmRNA[ele][0] > TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec][0] and  GffmRNA[ele][1] < TEEnds[rec][0] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "mRNA" + "_inside_TE"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("mRNA_inside_TE")
        elif GffmRNA[ele][0] < TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec][0] and  GffmRNA[ele][1] > TEEnds[rec][0] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "\t"+"TE_inside_"+ "mRNA" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_inside_mRNA")
        elif GffmRNA[ele][0] <= TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec][0] and  GffmRNA[ele][1] <= TEEnds[rec][0] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2Downstream_"+ "mRNA" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2Downstream_mRNA")
        elif GffmRNA[ele][0] >= TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec][0] and  GffmRNA[ele][1] >= TEEnds[rec][0] and TEStarts[rec][0] < GffmRNA[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2upstream_"+"mRNA" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2upstream_mRNA")
        elif GffmRNA[ele][0] == TEStarts[rec][0] and GffmRNA[ele][0] < TEEnds[rec][0] and  GffmRNA[ele][1] == TEEnds[rec][0] and TEStarts[rec][0] < GffmRNA[ele][1] :
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_exactly_"+ "mRNA"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_exactly_mRNA")
for ele in Gffexon:
    for rec in TEStarts:
        if Gffexon[ele][2]!= TEStarts[rec][1]:
            continue
        elif Gffexon[ele][0] > TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec][0] and  Gffexon[ele][1] < TEEnds[rec][0] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "exon" + "_inside_TE"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("exon_inside_TE")
        elif Gffexon[ele][0] < TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec][0] and  Gffexon[ele][1] > TEEnds[rec][0] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "\t"+"TE_inside_"+ "exon" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_inside_exon")
        elif Gffexon[ele][0] <= TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec][0] and  Gffexon[ele][1] <= TEEnds[rec][0] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2Downstream_"+ "exon" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2Downstream_exon")
        elif Gffexon[ele][0] >= TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec][0] and  Gffexon[ele][1] >= TEEnds[rec][0] and TEStarts[rec][0] < Gffexon[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2upstream_"+"exon" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2upstream_exon")
        elif Gffexon[ele][0] == TEStarts[rec][0] and Gffexon[ele][0] < TEEnds[rec][0] and  Gffexon[ele][1] == TEEnds[rec][0] and TEStarts[rec][0] < Gffexon[ele][1] :
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_exactly_"+ "exon"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_exactly_exon")
for ele in GffCDS:
    for rec in TEStarts:
        if GffCDS[ele][2]!= TEStarts[rec][1]:
            continue
        elif GffCDS[ele][0] > TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec][0] and  GffCDS[ele][1] < TEEnds[rec][0] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "CDS" + "_inside_TE"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("CDS_inside_TE")
        elif GffCDS[ele][0] < TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec][0] and  GffCDS[ele][1] > TEEnds[rec][0] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "\t"+"TE_inside_"+ "CDS" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_inside_CDS")	
        elif GffCDS[ele][0] <= TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec][0] and  GffCDS[ele][1] <= TEEnds[rec][0] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2Downstream_"+ "CDS" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2Downstream_CDS")
        elif GffCDS[ele][0] >= TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec][0] and  GffCDS[ele][1] >= TEEnds[rec][0] and TEStarts[rec][0] < GffCDS[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2upstream_"+"CDS" +"\t"+str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_fromInside2upstream_CDS")
        elif GffCDS[ele][0] == TEStarts[rec][0] and GffCDS[ele][0] < TEEnds[rec][0] and  GffCDS[ele][1] == TEEnds[rec][0] and TEStarts[rec][0] < GffCDS[ele][1] :
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_exactly_"+ "CDS"+"\t" +str(ele.rstrip("\n")))
            TEsinGENES[rec]=("TE_exactly_CDS")
for ele in Gfftranscript:
    for rec in TEStarts:
        if Gfftranscript[ele][2]!= TEStarts[rec][1]:
            continue
        elif Gfftranscript[ele][0] > TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec][0] and  Gfftranscript[ele][1] < TEEnds[rec][0] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "transcript" + "_inside_TE"+"\t" +str(ele.rstrip("\n"))+"\t"+rec.rstrip("\n").strip("\t")[3].strip("#")[1])
            TEsinGENES[rec]=("transcript_inside_TE")
        elif Gfftranscript[ele][0] < TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec][0] and  Gfftranscript[ele][1] > TEEnds[rec][0] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "\t"+"TE_inside_"+ "transcript" +"\t"+ ele)
            TEsinGENES[rec]=("TE_inside_transcript")
        elif Gfftranscript[ele][0] <= TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec][0] and  Gfftranscript[ele][1] <= TEEnds[rec][0] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2Downstream_"+ "transcript" +"\t"+ ele)
            TEsinGENES[rec]=("TE_fromInside2Downstream_transcript")
        elif Gfftranscript[ele][0] >= TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec][0] and  Gfftranscript[ele][1] >= TEEnds[rec][0] and TEStarts[rec][0] < Gfftranscript[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2upstream_"+"transcript" +"\t"+ ele)
            TEsinGENES[rec]=("TE_fromInside2upstream_transcript")
        elif Gfftranscript[ele][0] == TEStarts[rec][0] and Gfftranscript[ele][0] < TEEnds[rec][0] and  Gfftranscript[ele][1] == TEEnds[rec][0] and TEStarts[rec][0] < Gfftranscript[ele][1] :
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_exactly_"+ "transcript"+"\t" + ele)
            TEsinGENES[rec]=("TE_exactly_transcript")
for ele in GffcDNA_match:
    for rec in TEStarts:
        if GffcDNA_match[ele][2]!= TEStarts[rec][1]:
            continue
        elif GffcDNA_match[ele][0] > TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec][0] and  GffcDNA_match[ele][1] < TEEnds[rec][0] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "cDNA_match" + "_inside_TE"+"\t" +str(ele.rstrip("\n"))+"\t"+rec.rstrip("\n").strip("\t")[3].strip("#")[1])
            TEsinGENES[rec]=("cDNA_match_inside_TE")
        elif GffcDNA_match[ele][0] < TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec][0] and  GffcDNA_match[ele][1] > TEEnds[rec][0] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "\t"+"TE_inside_"+ "cDNA_match" +"\t"+ ele)
            TEsinGENES[rec]=("TE_inside_cDNA_match")
        elif GffcDNA_match[ele][0] <= TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec][0] and  GffcDNA_match[ele][1] <= TEEnds[rec][0] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2Downstream_"+ "cDNA_match" +"\t"+ ele)
            TEsinGENES[rec]=("TE_fromInside2Downstream_cDNA_match")
        elif GffcDNA_match[ele][0] >= TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec][0] and  GffcDNA_match[ele][1] >= TEEnds[rec][0] and TEStarts[rec][0] < GffcDNA_match[ele][1]:
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_fromInside2upstream_"+"cDNA_match" +"\t"+ ele)
            TEsinGENES[rec]=("TE_fromInside2upstream_cDNA_match")
        elif GffcDNA_match[ele][0] == TEStarts[rec][0] and GffcDNA_match[ele][0] < TEEnds[rec][0] and  GffcDNA_match[ele][1] == TEEnds[rec][0] and TEStarts[rec][0] < GffcDNA_match[ele][1] :
            print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) +"\t"+ "TE_exactly_"+ "cDNA_match"+"\t" + ele)
            TEsinGENES[rec]=("TE_exactly_cDNA_match")

for rec in TEStarts:
	if rec not in TEsinGENES:
		print (rec.rstrip("\n")+"\t"+str(int(TEsize[rec][1])-int(TEsize[rec][0])) + "TE_Not_in_genes"+ "-	--")

