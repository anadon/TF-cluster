########################################################################
@@@  @@@       @@@  @    @ @  @@@  @@@  @@@  @@@
 @   @         @    @    @ @  @     @   @    @ @
 @   @@@  @@@  @    @    @ @  @@@   @   @@@  @@
 @   @         @    @    @ @    @   @   @    @ @
 @   @         @@@  @@@  @@@  @@@   @   @@@  @ @
#TF-cluster#############################################################

##Table of Contents#####################################################

1............................................................Description

2..................................................................Usage

3...........................................................File Formats

4...............................................................Building

5.....................................................Build Requirements

6...........................................................Installation

7................................................................Credits

##Description###########################################################
TF-cluster is a tool made to cluster data, originally designed to find
interesting transcription factor co-expressions.

##Usage#################################################################
> ./triple-link-pthread -1 <FLOAT> -2 <FLOAT> -3 <FLOAT> -t <FILE PATH>
> -e <FILE PATH> -k <INTEGER> -c <"spearman" || "pearson"> 

Prints to stderr various status messages.  Results are printed to stdout
in the following format:
```
cluster 1:
<gene 1.1>
<gene 1.2>
...
<gene 1.M1>

cluster 2:
...

cluster N:
<gene N.1>
...
<gene N.MN>
```

##File Formats##########################################################

The settings file is no longer used, instead using command line 
arguments.

All input files expect UNIX standard line terminators.

The expression file used to generate the correlation matrix follows the
format below:
```
<NAME 1> <EXPRESSION COEFFICIENT 1> <EC 2> ... <EC M>
<NAME 2> <EXPRESSION COEFFICIENT 1> <EC 2> ... <EC M>
...
<NAME N> <EXPRESSION COEFFICIENT 1> <EC 2> ... <EC M>
```

The transcription file, used to specify transcriptions factors from the
expression file to create clusters from uses the following format:
```
<NAME 1>
<NAME 2>
...
<NAME N>
```

##Building##############################################################
> make

##Build Requirements####################################################
gcc-libs

libc

##Installation##########################################################
> make install

##Credits###############################################################
Programmer: Josh Marshall <jrmarsha@mtu.edu>

Professor:  Hairong Wei at hairong@mtu.edu
