########################################################################
@@@  @@@       @@@  @    @ @  @@@  @@@  @@@  @@@
 @   @         @    @    @ @  @     @   @    @ @
 @   @@@  @@@  @    @    @ @  @@@   @   @@@  @@
 @   @         @    @    @ @    @   @   @    @ @
 @   @         @@@  @@@  @@@  @@@   @   @@@  @ @
TF-cluster

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
> ./triple-link-pthread <PATH TO CONFIGURATION FILE>

Prints to stderr various status messages.  Results are printed to stdout
in the following format:

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


##File Formats##########################################################

The settings file acknowledges the following fields:

>expression=<PATH TO CORRELATION DATA>
>topPick=<NUMBER BETWEEN 1 and 255 inclusive>
>tripleLink1=<POSITIVE FLOATING POINT VALUE, LARGEST>
>tripleLink2=<POSITIVE FLOATING POINT VALUE, MIDDLE>
>tripleLink3=<POSITIVE FLOATING POINT VALUE, SMALLEST>

Settings which were acknowledged in the previous verion of TF-cluster,
but are no longer are:

>lib -- all resources are in the executable
>cpu -- the executable uses all available cores
>geneList -- this is no longer relevant
>kickSize -- The meaning of this has been opaque.  It may be re-added.

The expression file used to generate the correlation matrix follows the
format below:
><NAME 1> <EXPRESSION COEFFICIENT 1> <EC 2> ... <EC M>
><NAME 2> <EXPRESSION COEFFICIENT 1> <EC 2> ... <EC M>
>...
><NAME N> <EXPRESSION COEFFICIENT 1> <EC 2> ... <EC M>

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
