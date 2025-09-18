My output is as follows (located in my most recent commit under the Run Assembler Evaluation CL step).

Run ulimit -s 8192000
Dataset	Language	Runtime	N50
--------------------------------------------------------------------------------
data1	python	0:00:16	9990
data1	codon	0:00:10	9991
data2	python	0:00:34	9992
data2	codon	0:00:18	9993
data3	python	0:00:39	9824
data3	codon	0:00:18	9825
data4	python	0:07:43	159255
data4	codon	0:03:17	159478

As you can see, I made some modifications to the actions.yml file to farmilarize myself with the process and to help debug. I felt that the main, most important portion of this assignment was the familiarization.



BONUS QUESTION:

I used https://blast.ncbi.nlm.nih.gov/Blast.cgi, and fed it the contig.fasta files for data1, ... , data4. Here is what I found:

data1: A Porphyromonas gingivalis strain (specifically - best fit to "Porphyromonas gingivalis strain W83 chromosome, complete genome")

data2: Also a Porphyromonas gingivalis strain (best fit to "Porphyromonas gingivalis strain W83 chromosome, complete genome")

data3: A Paracidovorax citrulli strain (specifically - best fit to "Paracidovorax citrulli strain xjl12 chromosome, complete genome")

data4: I am getting the following error when trying to BLAST this data, so I was unable to decrypt this genome.
        Message ID#29 Error: Query string not found in the CGI context