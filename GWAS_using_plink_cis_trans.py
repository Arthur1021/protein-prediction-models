import argparse
import os, sys
import shutil
import multiprocessing.dummy as mp
parser = argparse.ArgumentParser()
parser.description='This script is used to run GWAS analysis with GLM for cis and trans region seperately.'
parser.add_argument('-in_gen',type = str,help=('input genotype file name with bfile (bed, bim, and fam) format, required'),required=True)
parser.add_argument('-TSS_file',type = str,help=('TSS file, required'),required=True)
parser.add_argument('-cis_range',type = int,default = 1000000,help=('define a cis region, based on TTS,default is 1000000'),required=False)
parser.add_argument('-pheno',type = str,help=('input file, phenotype residual, adjusted by covariate, required'),required=True)
parser.add_argument('-parallel',type = int,default = 30, help=('parallel computing, default is 30, optional'),required=False)
parser.add_argument('-pvalue',type = float,default = 5e-9, help=('p value filter, optional'),required=False)
parser.add_argument('-odir',type = str,default = './', help=('output dir'),required=False)
args = parser.parse_args()

# creat dir
if not os.path.exists(args.odir):
    os.mkdir(args.odir)
creat_dirs = []
creat_dirs.append(args.odir)
creat_dirs.append(args.odir+'/'+'phenotype_split')
creat_dirs.append(args.odir+'/'+'cis_geno')
# creat_dirs.append(args.odir+'/'+'trans_geno')
creat_dirs.append(args.odir+'/'+'cis_glm')
creat_dirs.append(args.odir+'/'+'trans_glm')
for item in creat_dirs:
    if not os.path.exists(item):
        os.mkdir(item)

#get FID
IID_to_FID ={}
infile = open (args.in_gen+'.fam','r')
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

#get cis region
range = args.cis_range

IDs = []
infile = open (args.TSS_file,'r')
head = infile.readline()
for myline in infile:
    myitem = myline.strip().split('\t')
    IDs.append(myitem[0])
    outfile = open (args.odir+'/cis_geno/'+myitem[0]+'_range.txt','w')
    if ',' not in myitem[1]:
        start = int(myitem[1].split(' ')[1])-range
        if start < 0:
            start = '0'
        else:
            start = str(start)
        end = str(int(myitem[1].split(' ')[1])+range)
        outfile.write(myitem[1].split(' ')[0]+'\t'+start+'\t'+end+'\t'+myitem[0]+'\n')
    elif ',' in myitem[1]:
        items = myitem[1].split(',')
        for item in items:
            start = int(item.split(' ')[1])-range
            if start < 0:
                start = '0'
            else:
                start = str(start)
            end = str(int(item.split(' ')[1])+range)
            outfile.write(item[0]+'\t'+start+'\t'+end+'\t'+myitem[0]+'\n')
     
    outfile.close()
infile.close()

def get_cis_region_SNP_genotype(ID):
    cmd = 'plink2 --bfile ' +args.in_gen+' --extract range '+args.odir+'/cis_geno/'+ID+'_range.txt --make-bed --out ' +args.odir +'/cis_geno/'+ID+'_cis'
    os.system(cmd)

def glm_cis(ID):
    cmd = 'plink2 --bfile '+args.odir+ '/cis_geno/'+ID+'_cis --pheno '+args.odir+ '/phenotype_split/'+ID+'.txt --glm hide-covar --out '+args.odir+'/cis_glm/'+ID +' --adjust '+"\n"
    os.system(cmd)

def glm_trans(ID):
    cmd = 'plink2 --bfile '+args.in_gen+ ' --pheno '+args.odir+ '/phenotype_split/'+ID+'.txt --glm hide-covar --out '+args.odir+'/trans_glm/'+ID + ' --pfilter '+str(args.pvalue)+"\n"
    os.system(cmd)

if __name__=="__main__":
    p=mp.Pool(args.parallel)
    p.map(get_cis_region_SNP_genotype,IDs)
    p.map(glm_cis,IDs)
    p.map(glm_trans,IDs)     
    p.close()
    p.join()
