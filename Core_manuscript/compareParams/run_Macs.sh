#/bin/bash

#$ -N macs2
#$ -P jknight.prjc -q short.qc
#$ -o stdout_Macs -e sterr_Macs.log -j y
#$ -cwd -V

### Run MACS2 On Set Of 3 BAM Files with User-Specified Parameters
##
## Idea is to take a directory containing 3 BAM files, representing
## technical replicates of an ATAC-seq experiment, and produce
## 3 corresponding NPF files. These files are named after the
## original BAM files, and deposited in an informatively
## named output directory...
##             ./macsq<$qVal>e<$extVal>_<sample>
## ... which contains the parameter values specified by the user
## for the MACS2 run. MACS2 sends it's outputs directly to this
## directory; run_Macs.sh can be run multiple times simultaneously 
## without issue. Two parameters are left to user control,the q-value, 
## which sets the FDR ($qVal) and an extsize, which controls 
## read length thereby influencing peak length ($extVal).

## In total, three inputs are required:
## $1 ---> directory containing BAM files
## $2 ---> q-value (cotrols the FDR)
## $3 ---> extsize (extension size, modulate read lengths)

## Since input BAM files lack a "chr" prefix for the chromosome 
## column, the last step of this script is to append "chr" to every
## row of the first column.

## Note on MACS2 parameters
## --shift : is set to be -1*half of extsize. This centers read density
## on the '5 cut points, and is recommended when calling peaks on
## DNase data, where the accessible region is defined by the cut point.
## --nomodel : because we are not doing ChIP, we do not need a model
## for shift size based on strand-specific peak-pairs.
## --nolambda : I feel this is an uncertain decision, but as we have
## no background sample (i.e. a ChIP input DNA), we choose not
## to use a local background calculated lambda.  In the absence of
## a background sample, MACS2 calculates the local lambda from the
## treatment itself, which means this value would be influenced by 
## actual signal, which is not in the spirit of a background threshold.

## JHendry, 2016/12/01


echo "**********************************************************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "**********************************************************************"


### Define Input Directory, Parameter Values,  Get Input File Name
bamDir=$1
bamFilePaths=$(ls $bamDir/*.bam)

qVal=$2
extVal=$3

echo "Input File Directory:"$bamDir
echo "No. Resident .bam Files:" `ls $bamDir | wc -w`
echo "File Names:"
printf "%s\n" $bamFilePaths
echo "Benjamini-Hochberg FDR (--qvalue):" $qVal
echo "Ext Size Value (--extsize):" $extVal
### Compute Shift Size From --extsize for DNase
# shift = -(1/2)*extsize
shiftVal=$[$extVal*-1/2]
echo "  Corresponding shift value (--shift):" $shiftVal 
echo ""

### Make Output Directory
npfDirPrefix=$(echo "macsq"${qVal#*.}"e"$extVal)
npfDirSuffix=${bamDir##*/}
npfDir=$(echo $npfDirPrefix"_"$npfDirSuffix)

if [ ! -d "$npfDir" ]; then
  mkdir $npfDir
  echo "Making Destination Directory:"$npfDir
  echo ""
fi

### Run MACS For All .bam Files in $bamDir
for bamFile in $bamFilePaths; do
  echo "Running MACS2 on " ${bamFile##*/}
    macs2 callpeak \
    --treatment $bamFile \
    --name ${bamFile##*/} \
    --outdir $npfDir \
    --nomodel \
    --nolambda \
    --extsize $extVal \
    --shift  $shiftVal \
    --qvalue $qVal \
    --verbose 3 
  echo "---------------------------------"
done

echo ""
echo "MACS2 Complete."
echo "Output Directory:" $npfDir
echo ""
echo "Appending 'chr' to chromosome designator (column 1)"
echo "in .narrowPeak files."

npfFiles=$(ls $npfDir/*.n*)
for npfFile in $npfFiles; do
  echo "Narrow Peak File Path:"
  echo $npfFile
  tempNpfFile=$(echo $npfFile".temp")
  echo "Temporary File Path"
  echo $tempNpfFile
  cat $npfFile | awk '{$1="chr"$1; print}' > $tempNpfFile && mv $tempNpfFile $npfFile 
  echo "Append complete."
done



echo "**********************************************************************"
echo "Finished at: "`date`
echo "**********************************************************************"
echo ""
echo ""
echo ""
