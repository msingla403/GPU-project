# GPU-project
Course project

#### Compile Instructions
```console
$ nvcc main.cu
```

#### Run 
```consile
$ ./a.out dataset(sparse).txt 
```

### Instructions
Sparse matrix input format is as follows:
1. The first line contains 3 integers, num-rows, num-columns, num-entries
2. The following num-entries lines contains 2 integers, row-id and column-id

These are 1-indexed, also we do not need value at that position in matrix.

There are 5 parmaeters that can be tuned to experiment for getting better results
1. PANEL SIZE
2. DENSE THRESHOLD
3. SIGNATURE LENGTH
4. BAND SIZE
5. NUMBER OF HASH BUCKETS

There is also a DEBUG flag to print helpful arrays computed in runtime.
