## Multi-Array Correl. Computation Utility
This respository is currently for describing the `OneStopWrapper.pl` of [maccu](https://maccu.sourceforge.net/index.html), which performs co-expression clustering in one command. In this document, we also present the use of a co-expression database of Arabidopsis Col-0 RNAseq samples made based on data from [DEE2](https://dee2.io/).

### Installation requirement
It should be feasible to run `OneStopWrapper.pl` in recent linux distributions. The only requirement is java runtime environment greater than or equal to 8. In this document, we use an ubuntu20 VM with 8GB ram to run all commands.

Install java runtime environment.
```
ubuntu@maccu:~$ sudo apt update
ubuntu@maccu:~$ sudo apt install openjdk-11-jre-headless
```

Downlaod the wrapper package and extract it. Must be sure to keep extracted files been placed under the same directory.
```
ubuntu@maccu:~$ wget https://downloads.sourceforge.net/project/maccu/0.8/maccuWrapper.tar.gz
ubuntu@maccu:~$ tar -zxvf maccuWrapper.tar.gz
```

### Download the co-expression database
Here we provide the co-expression database made based on data from the (DEE2)[https://dee2.io/]. Please note that this database and the DEE2 database are both under [GPLv3 public license](https://www.gnu.org/licenses/gpl-3.0.en.html) and allow academic and commercial use. We thank DEE2 and ask you to consider supporting DEE2 if your organization uses DEE2 data for financial benefit.

Download the co-expression database and extract it.
```
ubuntu@maccu:~$ wget https://maccu.project.sinica.edu.tw/maccu/coexDB/coexDB20230714.tar.gz
ubuntu@maccu:~$ tar -zxvf coexDB20230714.tar.gz
coexDB20230714/
coexDB20230714/ath/
coexDB20230714/ath/Col0/
coexDB20230714/ath/Col0/README.txt
coexDB20230714/ath/Col0/sel20210116.col0.TMM.ALL
coexDB20230714/ath/Col0/sel20210116.col0.TMM.flower
coexDB20230714/ath/Col0/sel20210116.col0.TMM.leaf
coexDB20230714/ath/Col0/sel20210116.col0.TMM.root
coexDB20230714/ath/Col0/sel20210116.col0.TMM.rosette
coexDB20230714/ath/Col0/sel20210116.col0.TMM.seed
coexDB20230714/ath/Col0/sel20210116.col0.TMM.seedling
coexDB20230714/ath/Col0/sel20210116.col0.TMM.shoot
coexDB20230714/ath/Col0/sel20210116.col0.TMM.whole
```

The files named `sel20210116.col0.TMM.*` are our database files in tab-delimited text format. They are all read count matrix normalized by the TMM method (PMID: 20196867). Their columns are samples and rows are genes. The one suffixed by `ALL` are composed of 5556 Arabidopsis Col-0 RNAseq samples, which were selected from the DEE2 database following a series of considerations. All other database files are extracted portions of this `ALL` file, where the extracted portions were decided by parsing metadata download from NCBI. For example, the `root` one should be composed of root-related samples. For another example, the `whole` one should be composed of samples using _whole plants_.

### First execution of `OneStopWrapper.pl`
It is OK to execute the script by specifying the path, and it is also OK to put its path into the PATH environment variable. Simple description of options will be displayed if no options entered.
```
ubuntu@maccu:~$ maccuWrapper/OneStopWrapper.pl
Usage: OneStopWrapper.pl [<options>]+ <tgzFile>
 options:
   [-input <name> <fileName>]+ : input gene set names and gene set files, at least one.
   -inputParam <parameterStr>  : fishing/clustering option (default: -L true)
   -maccu <JARpath>            : path to maccu.jar (default: maccuWrapper/maccu.jar)
   -method <classStr>          : maccu.CoExpressFishing or maccu.RelationComputer(default)
   [-dbase <DBname> <DBloc>]+  : coex database names and database files
   -probeMap <fileName>        : probe-gene mapping file, in case of using array database (default: /dev/null)
   -thresh <from> <to> <step>  : correlation threshold range and step (default: 0.70 0.95 0.05)
   -outDir <dirName>           : output folder name in the final tar.gz file (default: coex)
   -graphAdj <optionStr>       : post computation adjustment parameters (default: "")
   [-fold <name> <foldChangeFileName>]+ : fold-change set name and files (default: no setting)
   [-graphAdjStr <inputPrefix> -remove <operandPrefix> <outputPrefix>]+ : graph-level remove operation on graphs
                               made base on different coex databases
   -graphAdjShift <shift>      : graph-level operand shift on computed network arrays (default: 1)
   -keep                       : keep working directory (default: no keep)
   -outputTGZ <outputFilename> : output .tgz filename (default: coexRes.tgz)
```

### Download example data files and run the first example
This github repository contains example files under the `example` folder.
```
ubuntu@maccu:~$ git clone https://github.com/wdlingit/maccu.git

ubuntu@maccu:~$ cd maccu/example/deg/
ubuntu@maccu:~/maccu/example/deg$ ls
foldchange.txt  genelist.txt
```

The `deg` directory contains two files, they are from experiment data of PMID: 20181752. The `genelist.txt` file contains simply gene IDs in AGI accessions. The `foldchange.txt` file is in tab-delimited format, where the first column for gene IDs and the second column for fold-chagne values. It is strongly suggested to use log2-fold-change as the fold-change values because that could be handled easily in Cytoscape.

Since the experiment is for studying Arabidopsis roots and the initial analysis (of array) discovered 187 genes, it is of our interests to figure out potential co-expression modules inside the 187 genes. To do that, we use the database of root-related samples to infer potential co-expression modules.
```
ubuntu@maccu:~/maccu/example/deg$ /home/ubuntu/maccuWrapper/OneStopWrapper.pl \
                                  -input DEG genelist.txt \
                                  -dbase root /home/ubuntu/coexDB20230714/ath/Col0/sel20210116.col0.TMM.root \
                                  -fold foldchange foldchange.txt
```
