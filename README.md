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
