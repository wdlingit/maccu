## Multi-Array Correl. Computation Utility
This respository is currently for describing the `OneStopWrapper.pl` of [maccu](https://maccu.sourceforge.net/index.html), which performs co-expression clustering in one command. In this document, we also present the use of a co-expression database of Arabidopsis Col-0 RNAseq samples made based on data from [DEE2](https://dee2.io/).

### Installation requirement
It should be feasible to run `OneStopWrapper.pl` in recent linux distributions. The only requirement is java runtime environment greater than or equal to 8. In this document, we use an ubuntu20 VM with 8GB ram for running all commands.

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

The files named `sel20210116.col0.TMM.*` are our database files in tab-delimited text format. They are all read count matrix whose columns are samples and rows are genes. The one suffixed by `ALL` are composed of 5556 Arabidopsis Col-0 RNAseq samples, where were selected from the DEE2 database following a series of considerations. All other database files are extracted portions of this `ALL` made by parsing metadata download from the NCBI database. For example, the `root` one should be composed of root-related samples. For another example, the `whole` one should be composed of samples using _whole plants_.
