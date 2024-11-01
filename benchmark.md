```shell
hyperfine --warmup 3 --parameter-list threads 2,4,6,8,16,32,64 './target/release/star-aligner {threads} > /dev/null' '~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN {threads} --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null' --export-markdown benchmark.md
```

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `./target/release/star-aligner 2 > /dev/null` | 1.122 ± 0.015 | 1.104 | 1.150 | 2.48 ± 0.10 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 2 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 1.991 ± 0.005 | 1.985 | 2.000 | 4.40 ± 0.16 |
| `./target/release/star-aligner 4 > /dev/null` | 0.637 ± 0.003 | 0.633 | 0.642 | 1.41 ± 0.05 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 4 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 2.080 ± 0.004 | 2.075 | 2.088 | 4.59 ± 0.17 |
| `./target/release/star-aligner 6 > /dev/null` | 0.505 ± 0.010 | 0.496 | 0.525 | 1.11 ± 0.05 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 6 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 2.172 ± 0.005 | 2.164 | 2.179 | 4.80 ± 0.18 |
| `./target/release/star-aligner 8 > /dev/null` | 0.453 ± 0.017 | 0.437 | 0.499 | 1.00 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 2.256 ± 0.003 | 2.251 | 2.260 | 4.98 ± 0.18 |
| `./target/release/star-aligner 16 > /dev/null` | 0.453 ± 0.026 | 0.412 | 0.482 | 1.00 ± 0.07 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 2.605 ± 0.003 | 2.600 | 2.609 | 5.75 ± 0.21 |
| `./target/release/star-aligner 32 > /dev/null` | 0.606 ± 0.020 | 0.563 | 0.631 | 1.34 ± 0.07 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 32 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 3.350 ± 0.004 | 3.344 | 3.358 | 7.40 ± 0.27 |
| `./target/release/star-aligner 64 > /dev/null` | 1.479 ± 0.235 | 1.095 | 1.835 | 3.27 ± 0.53 |
| `~/STAR-2.7.2a/bin/Linux_x86_64/STAR --runThreadN 64 --genomeDir /data/wenjie/star-aligner/star-aligner/test/ercc92-1.2.0/star/ --readFilesIn ~/test.2w.fastq  --sjdbOverhang 100 --outFilterScoreMin 0 --readNameSeparator space --outSAMunmapped Within KeepPairs --outSAMtype SAM  --outStd SAM --outSAMorder PairedKeepInputOrder > /dev/null` | 4.964 ± 0.021 | 4.945 | 5.005 | 10.96 ± 0.41 |
