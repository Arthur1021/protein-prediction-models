import argparse
import glob
import os, sys
import shutil
import multiprocessing.dummy as mp
parser = argparse.ArgumentParser()
parser.description='This script is used to build model and perform fusion analysis.'
parser.add_argument('-in_gen',type = str,help=('input genotype file name with bfile (bed, bim, and fam) format,used for generating 1000g LD list file, required'),required=True)
parser.add_argument('-pheno',type = str,help=('input file, phenotype file, required'),required=True)
parser.add_argument('-bfdir',type = str,help=('input dir, extract up and down 100k bfile, used for build model, required'),required=True)
parser.add_argument('-parallel',type = str,default = 30, help=('parallel computing, default is 30, optional'),required=False)
parser.add_argument('-odir',type = str,default = './', help=('output dir, default is ./, do not change it'),required=False)
args = parser.parse_args()

#creat dirs
if not os.path.exists(args.odir):
    os.mkdir(args.odir)

pheno_dir = args.odir+'/phenotype_split/'
if not os.path.exists(pheno_dir):
    os.mkdir(pheno_dir)

os.system("mkdir -p phenotype_split")
os.system("mkdir -p TWAS_fusion/tmp/")
os.system("mkdir -p TWAS_fusion/out/")
os.system("mkdir -p output/TWAS_fusion/tmp/")
os.system("mkdir -p output/TWAS_fusion/out")


#if break, continue
exist_result = {}
exist_results = glob.glob(args.odir+'/TWAS_fusion/out/*.100kb.wgt.RDat')
for myfile in exist_results:
    id = myfile.split('/')[-1].replace('.all.100kb.wgt.RDat','')
    exist_result[id] = 1

path = args.bfdir
files = os.listdir( path )
full_path = os.path.abspath(path)
fbile = []
for myfile in files:
    if "unambig.Z.bed" in myfile and myfile.replace('.noDup.unambig.Z.bed','') not in exist_result:
        fbile.append(myfile[:-4])
fbile.sort()

IID_to_FID = {}
infile = open (args.bfdir+'/'+fbile[0]+'.fam','r')
for myline in infile:
    myitem = myline.strip().split()
    IID_to_FID[myitem[1]] = myitem[0]
infile.close()

#split phenotype
def split_pheno_file(pheno_file):
    infile = open (pheno_file,'r')
    phenotype_split ={}
    head = infile.readline().strip().split('\t')
    for myline in infile:
        myitem = myline.strip().split('\t')
        for i in range(1,len(myitem)):
            phenotype_split.setdefault(head[i],[]).append(IID_to_FID[myitem[0]]+'\t'+myitem[0]+'\t'+myitem[i])
    infile.close()
    return phenotype_split
result = split_pheno_file(args.pheno)

for k,v in result.items():
    outfile = open (args.odir+'/phenotype_split/'+k+'.txt','w')
    outfile.write('FID'+'\t'+'IID'+'\t'+k+'\n')
    for item in v:
        outfile.write(item+'\n')
    outfile.close()



def model_TWAS_fusion(myfile):
    name = myfile.split("/")[-1].replace('.noDup.unambig.Z','')
    print ("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.compute_weights.modified.R --bfile " +full_path+"/"+ myfile + " --tmp TWAS_fusion/tmp/" + name+ ".all.100kb --out TWAS_fusion/out/" + name +".all.100kb --verbose 2 --pheno phenotype_split/" + name+ ".txt --noclean --PATH_gcta /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/gcta_nr_robust --PATH_gemma /mnt/lvm_vol_1/dghoneim/EGA/gemma-0.98.1-linux-static")
    os.system("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/FUSION.compute_weights.modified.R --bfile " +full_path+"/"+ myfile + " --tmp TWAS_fusion/tmp/" + name+ ".all.100kb --out TWAS_fusion/out/" + name +".all.100kb --verbose 2 --pheno phenotype_split/" + name+ ".txt --noclean --PATH_gcta /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/gcta_nr_robust --PATH_gemma /mnt/lvm_vol_1/dghoneim/EGA/gemma-0.98.1-linux-static")
    print ("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/read.wgt.RDat.R TWAS_fusion/out/"+ name+".all.100kb.wgt.RDat TWAS_fusion/out/" + name +".all.100kb.wgt.RDat.ld.list")
    os.system("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/read.wgt.RDat.R TWAS_fusion/out/"+ name+".all.100kb.wgt.RDat TWAS_fusion/out/" + name +".all.100kb.wgt.RDat.ld.list")
    print ("plink --bfile " + args.in_gen + " --allow-extra-chr --extract TWAS_fusion/out/"+name+ ".all.100kb.wgt.RDat.ld.list --make-bed --out " +name +".ld.1000g.Z")
    os.system("plink --bfile " + args.in_gen + " --allow-extra-chr --extract TWAS_fusion/out/"+name+ ".all.100kb.wgt.RDat.ld.list --make-bed --out " +name +".ld.1000g.Z")
    print ("echo \"WGT     ID      CHR     P0      P1\" > " +name+".pos")
    os.system("echo \"WGT     ID      CHR     P0      P1\" > " +name+".pos")
    print ("echo \""+name+".all.100kb.wgt.RDat       "+name+"    Z       1       9999999999\" >> " +name +".pos")
    os.system("echo \""+name+".all.100kb.wgt.RDat       "+name+"    Z       1       9999999999\" >> " +name +".pos")



if __name__=="__main__":
    p=mp.Pool(int(args.parallel))
    p.map(model_TWAS_fusion,fbile)
    p.close()
    p.join()
