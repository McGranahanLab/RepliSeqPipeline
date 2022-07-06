# RepliSeqPipeline
Pipeline to process fastq files provided by the Repli-Seq protocol of G1, early S-phase and late S-phase reads.

# Getting started & Documentation
1) Before running the pipeline a reference fasta file (e.g. ucsc.hg19.fasta) and a corresponding file with chromsome sizes (hg19.chr.sizes.txt) needs to be downloaded.

2) The pipeline requires a design file with the following columns as input (see examples/designFile.txt):

      ```
      Patient = patient id
      Sample = combined name of Patient and Type
      Type = FACS sorted sample (G1, Early or Late)
      FQ1 = path to fastq file
      FQ2 = path to fastq file if paired-end sequencing was performed, NA otherwise
      ```

3) Adapt parameters and paths in RepliSeq.wrapper.sh file

More details about the pipeline itself can be found in the pipelineOverview.pptx file and the corresponding publication xxx.
