# Assignment the First

## Part 1
1. Be sure to upload your Python script (i, ii, and iii).

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
![](https://github.com/christian-lafrance/Demultiplex/blob/master/Assignment-the-first/Index_1.png)
![](https://github.com/christian-lafrance/Demultiplex/blob/master/Assignment-the-first/Index_2.png)
![](https://github.com/christian-lafrance/Demultiplex/blob/master/Assignment-the-first/Read_1.png)
![](https://github.com/christian-lafrance/Demultiplex/blob/master/Assignment-the-first/Read_2.png)


    2. 30 seems like a good cut off. All of the average quality scores
        are above 30. A score of 30 would indicate a 0.1% chance of an
        incorrect base call. 
    3. 
        Index_1 file:
        '''zcat 1294_S1_L008_R2_001.fastq.gz | sed -n "2~4p" | grep "N" | wc -l'''
        Result: 3976613 reads contain an N. 
        Index_2 file:
        '''zcat 1294_S1_L008_R3_001.fastq.gz | sed -n "2~4p" | grep "N" | wc -l'''
        Result: 3328051 reads contain an N. 
    
## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
