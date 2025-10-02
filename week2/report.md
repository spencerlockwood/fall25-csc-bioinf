This (Week 2) deliverable focused on code analysis and motif parsing, with the goal being porting a motif discovery pipeline from a Python package (BioPython or TRviz) into Codon, and creating tests that can run both in Python and Codon environments.

I ran into many difficulties using Biopython, so I decided to use TRviz, which turned out successful.

The first step of my process was to convert decomposer.py and utils.py (from TRviz) from python to Codon. I made sure to skip visualizer.py and other non-essential parts, here.

Then, I created a motifs.py file based on that of which is in Biopython, so I may test accurately.

Then, I created test.py with the necessary tests to prove my importing was successful. From here I created my actions.yml file and pushed my finished product.

My main struggle in this deliverable was creating the motifs.py file that worked successfully.