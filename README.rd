<<RBSeval>>

 RBSeval predicts protein abundance from mRNA sequence around the start codon.
 In this method, protain abundance is estimated based on the accessibility around
 the start coodon and the predicted activity of the Shine-Dalgarno sequence.  
 The former and latter is calculated by the Raccess program [1] and the EMOPEC method [2],
 respectively. RBSeval combines the above two features using the linear regression and 
 outputs a predicted protein abundance value.

<<Prerequisite>>
 RBSeval uses the Raccess program for calculating accessibility around the start
 codon. You must have the Raccess program enabled before using the RBSeval.
 Raccess is currently available from "https://github.com/gterai/raccess".

 RBSeval uses the EMOPEC method to estimate the strength of the Shine-Dalgarno (SD) site.
 Before using RBSeval, you need to have the EMOPEC Python library installed. EMOPEC is
 available at “https://github.com/micked/EMOPEC”. You don’t need all of EMOPEC’s dependenci
 the following command is enough to use RBSeval:
 pip install git+https://github.com/micked/EMOPEC.git

<<Usage>>
 perl RBSeval.pl [5'-UTR (fasta file)] [CDS (fasta file)] [run_raccess_contrafold (binary file)]

<<Example>>
 perl RBSeval.pl example/test_utr.fasta  example/test_cds.fasta /your_path_to_raccess/src/raccess/run_raccess_contrafold

<<Output>>
 accC            : -2.67845
 EMOPEC score    : 1.62256087658972
 predicted SD    : ACGACA
 RBSeval score   : 7.16090606082678
 Exp. level      : Very low (Bottom 30 percentile)

 "accC" is the accessibility around the start codon. "EMOPEC sore" is the predicted activity of 
 the Shine-Dalgarno sequence. "predicted SD" is the plausible SD sequence predicted by the
 EMOPEC method. "RBSeval score" is a predicted protein abundance value. "Exp. level" is the
 protein expression level assigned based on the RBSeval scores for a dataset used to train
 RBSeval (see [3] for details of the dataset). 

<<Citation>>
 Please site [3] when you use RNAeval.

<<Reference>>
[1] A detailed investigation of accessibilities around target sites of siRNAs and miRNAs.
    Kiryu H, Terai G, Imamura O, Yoneyama H, Suzuki K, Asai K.
    Bioinformatics. 2011 Jul 1;27(13):1788-97.
[2] Predictable tuning of protein expression in bacteria.
    Bonde MT, Pedersen M, Klausen MS, Jensen SI, Wulff T, Harrison S, Nielsen AT, Herrgård MJ,
    Sommer MO.
    Nat Methods. 2016 Mar;13(3):233-6.
[3] Improving the prediction accuracy of protein abundance in Escherichia coli using mRNA accessibility (submitted).
