# splitter tools

Building splitter tools works via cmake. Therefore please perform the following steps:

1. `mkdir build`
2. `cd build`
3. `cmake ..` - make sure that cmake sees at least gcc/g++ 4.7.4
4. `make & make install`

After `make install` you will have an include, lib and bin folder in your root directory.

Within the bin folder you can use the tools.

`gtxtools -gtx gtf/gtt-file -stats statsfile`

where you stats-file is a tab-seperated file on the features you want statistics for:

chromosome	gene	0
gene	transcript	0
transcript	exon	0
transcript	exon	1
transcript	exon	2

There are 3 different modes:

0. calculates statistics on any gene child feature of chromosome
1. calculates statistics on any non-gene child interval of chromosome
2. as 1. but does not include first/last interval
