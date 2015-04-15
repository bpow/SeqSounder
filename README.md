# SeqSounder
A 'sounder' to determine depth of coverage of genetic sequence. This is really just a test, but may be useful to
someone.

# Building
1. Checkout from github
2. `./gradlew jar` from within the main directory. This will produce a SeqSounder-VERSION.jar file

# Using
Try `java -jar SeqSounder-VERSION.jar -h` for a list of options. Outputs will be a 'report' file, and, optionally,
a BedGraph-formatted file or 'coverage fasta' file. You have to provide at least one bam file, of course...

If you don't specify any intervals (with `-l` or `-r`), then the output files will get quite large. Trade disk space
for CPU with the `-z` option to compress the larger files as they are written.

# Extracting regions below a certain coverage
I prefer to use the "BedGraph" output for this. For a bamfile named `sample.bam`, with target intervals of `intervals.bed`, if you run:

`java -Xmx2g -jar SeqSounder-VERSION.jar -l intervals.bed -s suffix -t sample.bam`

... then you will get files named `sample.suffix.bedgraph` and `sample.suffix.report` in the same directory as your bam file (or a `sample.suffix.bedgraph.gz` if you had used the `-z` option.

Extracting regions with a coverage of less than 20 is as simple as:

`awk '$4 < 20' <sample.suffix.bedgraph`

If you want the total size of these regions, you can use the more complicated awk expression:

`awk '$4 < 20 {SUM+=$3-$2} END {print SUM}' <sample.suffix.bedgraph`
