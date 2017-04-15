# phyNGSC

phyNGSC is a hybrid strategy between MPI and OpenMP to accelerate the compression of big FASTQ datasets by combining the best features of distributed and shared memory architectures to:

* Balance the load of work among processes.
* Alleviate memory latency by exploiting locality.
* Accelerate I/O by reducing excessive read/write operations and inter-node message exchange. 

Our algorithm introduces a novel timestamp-based approach which allows concurrent writing of compressed data in a non-deterministic order and thereby allows us to exploit a high amount of parallelism.

As a proof-of-concept, we implemented some methods developed for [DSRC v1](http://sun.aei.polsl.pl/dsrc/) to underline the compression portion of our hybrid parallel strategy, since it exhibits superior performance for sequential solutions.

### Disclaimer:

All files provided in this project are part of the article "*[A Hybrid MPI-OpenMP Strategy to Speedup the Compression of Big Next-Generation Sequencing Datasets](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7895161&isnumber=4359390)*" accepted for publication in the **IEEE Transactions on Parallel and Distributed Systems** journal and should be used for testing purposes.

## Getting Started

### Prerequisites

In order to compile, build and run phyNGSC you need to install and configure the following:

* __OS:__ Any Linux distribution, preferably Ubuntu (14.04 LTS or a later) or CentOS.
* __Compiler:__ The GNU Compiler Collection (GCC) version 4.8 or later.
* __MPI__: MPICH version 3 or later. MVAPICH2, the Ohio State University derivative of MPICH, can also be use.

### Compiling and Building

To compile an build phyNGSC we provided a simple `Makefile`. Navigate to the directory containing the files and run the command:

```
make
```

If your system configurations (compiler and MPI implementation) are correct, you should see the following output:

```
mpicxx -O3 -m64  -Wall -fopenmp -std=c++11 -c phyNGSC.cpp -o phyNGSC.o
mpicxx -O3 -m64  -Wall -fopenmp -std=c++11 -c bit_stream.cpp -o bit_stream.o
mpicxx -O3 -m64  -Wall -fopenmp -std=c++11 -c huffman.cpp -o huffman.o
mpicxx -O3 -m64  -Wall -fopenmp -std=c++11 -c tasks.cpp -o tasks.o
mpicxx -O3 -m64  -Wall -fopenmp -std=c++11 -o phyNGSC phyNGSC.o bit_stream.o huffman.o tasks.o
```

To clean up object files generated during compilation in the directory and the `phyNGSC` program, run:

```
make clean
```

## Running phyNGSC

To run phyNGSC, execute the command

```
mpiexec -np <num_of_processes> ./phyNGSC <input_file> <output_file> <num_of_threads>
```
where

* `mpiexec` – Runs an MPI program.
* `-np <num_of_processes>` – Specify the number of processes to use. `<num_of_processes>` must be at least 2.
* `./phyNGSC` – Executable phyNGSC program.
* `<input_file>` – A FASTQ file with file extension `*.fastq`. Must exist.
* `<output_file>` – An NGSC file with file extension `*.ngsc`. Created if it does not exist.
* `<num_of_threads>` – Number of OpenMP threads to use per MPI process. `<num_of_threads>` must be at least 1.

**NOTE:** Depending on the implementations of MPI currently on your system, `mpicxx` (to compile and build) and `mpiexec` (to run) can be wrappers for MPICH or OpenMPI (to mention two MPI implementations). Make sure `mpicxx` and `mpiexec` are using MPICH.

### Usage Example:
Compress the file `input10MB.fastq` using 8 processes and 5 thread per process, and name the resulting NGSC file `output.ngsc`:

```
mpiexec -np 8 ./phyNGSC input10MB.fastq output.ngsc 5
```

The output should look similar to the one below:

```
[I] INFO: Processing <<input10MB.fastq>> with 8 MPI processes and 5 threads per process.

RANK   COMP_TIME   N_BLOCK   N_SUBBLOCKS
----------------------------------------------
0      3.269169    3         16
2      3.101142    3         16
3      3.101125    2         16
4      3.268828    2         16
5      3.268804    3         16
6      3.211773    3         16
1      3.143939    3         16
7      3.211767    3         16
```

## Dataset Files

In this project a FASTQ file named `input10MB.fastq` is provided for your testing purposes. Each record in a FASTQ file is compose by 4 lines and phyNGSC can only recognize FASTQ files with the following record structure:

1. Title line starting with the **@** sign.
2. DNA sequence line.
3. Third line of the record starting with **+** sing and **NO** repetition of the title line.
4. Quality score line should be the same length of the corresponding DNA sequence line. 

If the FASTQ file contains records with the title line repeated in the third line, phyNGSC will produce an error.

### Example of Valid FASTQ Files for phyNGSC
The first 5 records of the file `input10MB.fastq` are shown bellow:
```
@ERR005195.1 BGI-FC30BFTAAXX_5_1_000:1689/2
CTCCCATATCCTTAGAGAAAATCCCCAATGCCTAGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIII'8=;I?DG&
@ERR005195.2 BGI-FC30BFTAAXX_5_1_000:125/2
TGTTTGGCAAGGTCCTACAAAAGTTGCAACTCTCAC
+
IIIIIIIIIIIIIIIIIIII2?9IIII*7)IIII'-
@ERR005195.3 BGI-FC30BFTAAXX_5_1_000:137/2
CAAGGCAGGCGGGTCACTTGAGGTCAGGAGTTTGAG
+
IIIIIIIIIIIIII@:IIII+9I*7+5I(/22I2#>
@ERR005195.4 BGI-FC30BFTAAXX_5_1_000:1108/2
TTACTCACGACCAAGCGTCATTTTCCATTATTGCTG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII&II(
@ERR005195.5 BGI-FC30BFTAAXX_5_1_000:477/2
TTCTAAGCCACTGGGGAATTGATGGATGCTCCCAGG
+
IIIIIIIIIIIIIIIIIIIIIIIE6III9IIII4:(
```

### Example of Invalid FASTQ Files for phyNGSC
These are the first 5 records of a FASTQ file containing repeated title in the third line of each record:

```
@SRR013667.1 30PTAAAXX:5:1:0:1203 length=76
NCCAGCAGCCATAACTGGAATGGGAAATAAACACTATGTTCAAAGCAGAGAAAATAGGAGTGTGCAATAGACTTAT
+SRR013667.1 30PTAAAXX:5:1:0:1203 length=76
#>A@BABAAAAADDEGCEFDHDEDBCFDBCDBCBDCEACB>AC@CDB@>>CB?>BA:D?9>8AB685C26091:77
@SRR013667.2 30PTAAAXX:5:1:0:502 length=76
NGGAAAGAAATGAAATGGGATGGAACAACCCGAAGGGAAGGGAAGGAAATGGAGAGTAAGGGAGCTGACCAGTATC
+SRR013667.2 30PTAAAXX:5:1:0:502 length=76
#>A@B?BAA&AADDEGCEFDHDE.BCF)B)DB;BDCEAC+>A+@BDD@C>CB?>:*:D.9>8=B.A5;+/09*477
@SRR013667.3 30PTAAAXX:5:1:0:703 length=76
NACAGGAAGAGCGATTCACTTAGGATGACAGGAACTCAAATCTACGCTTGATGCAGACTATACTATGTGATATCAA
+SRR013667.3 30PTAAAXX:5:1:0:703 length=76
#>A@BABAAAAADAEGCEF9HDE+BC.D7CE(CB)*CACB>'CDC0B>,1(B?>:.:DB9;83B6?6;26/<149;
@SRR013667.4 30PTAAAXX:5:1:0:407 length=76
NAGCTGCCATTTTCTCATCTGTGAACTGGAATGATAGATACCACTCCACATCCCATCAAGGGTCAATGAACAGTAG
+SRR013667.4 30PTAAAXX:5:1:0:407 length=76
#>A@BABAAAAADDEGCEFDHDEDBCFDBCDBCBDCEACB>AC@CDB@>>CB?>:A:D?9B8<B685;26091472
@SRR013667.5 30PTAAAXX:5:1:0:1139 length=76
NGAAAACAGCAACAAACGACCGAGACAGACGACCACTAGTACTGACACGAATGTTGAGACTTCAAAGGACAAAACC
+SRR013667.5 30PTAAAXX:5:1:0:1139 length=76
#>A@BABAAAADDDEGC7FD+DEBBCFDB1DB<EDC(AC9>AC@CDB@>=C(?.:A:DC9E8=(38C:2:B918?5
```

### More Datasets for Testing

* [ERR229788.filt.fastq](https://goo.gl/6y1Qw3)
* [ERR260401_2.filt.fastq](https://goo.gl/oYV2Zq)
* [ERR022729_1.filt.fastq](https://goo.gl/QoDxen)
* [ERR022548_2.fastq](https://goo.gl/Zqw9WT)

You can find more FASTQ datasets at the [EBI FTP site](https://goo.gl/CkWqEY) or the [IGSR data portal](http://www.internationalgenome.org/data-portal/sample).

## Authors

* **Sandino Vargas-P&eacute;rez** - Department of Computer Science at Western Michigan University, USA. [E-mail](mailto:sandinonarciso.vargasperez@wmich.edu).
* **Fahad Saeed** - Department of Electrical and Computer Engineering and Department of Computer Science, Western Michigan University, USA. [E-mail](mailto:fahad.saeed@wmich.edu).

## Acknowledgment

This work was supported in part by grant NSF CRII CCF-1464268. The authors would like to acknowledge [XSEDE](https://www.xsede.org) startup grant CCR150017 and the help of staff at the SDSC Gordon computer cluster. We would also want to acknowledge Sebastian Deorowicz and Szymon Grabowskie for their DSRC algorithm (version 1) which underlines the compression portion of our parallel hybrid strategy.