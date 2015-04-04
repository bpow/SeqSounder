package cser;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPOutputStream;

public class CallableDepth {
    private static final int COVERAGE_HISTOGRAM_MAX = 1000;
    @Parameter(names = {"-h", "--help"}, help = true)
    public boolean help;
    @Parameter(names="-mapQ", description="Minimum mapping quality")
    public int minimumMapScore = 20;
    @Parameter(names="-baseQ", description="Minimum base quality")
    public int minimumBaseQuality = 20;
    @Parameter(names="-keepDupes", description="include duplicates in analysis")
    public boolean keepDupes;
    @Parameter(names={"-z", "-compress"}, description="Compress output files while writing")
    public boolean compressOutput;
    @Parameter(names={"-bed", "-l"}, description="BED file of intervals to analyze (ZERO-based coordinates)")
    public String bedFile = null;
    @Parameter(names={"-r", "-region"}, description="Regions to analyze (ONE-based), like \"1:32221-42212\", can be specified multiple times")
    public List<String> region = new ArrayList<String>();
    @Parameter(names={"-f", "-covFasta"}, description="Generate a .covfasta file")
    public boolean makeCovFasta;
    @Parameter(names={"-t", "-covBedGraph"}, description="Generate a tab-delimited (bedGraph) coverage file")
    public boolean makeCovBedGraph;
    @Parameter(names={"-s", "-suffix"}, description="Additional suffix to add to output files")
    public String suffix = "";
    @Parameter(description = "bamFiles")
    public List<String> bamFiles = new ArrayList<String>();

    private ArrayList<Interval> setupIntervals() throws IOException {
        ArrayList<Interval> inputIntervals;
        if (bedFile != null) {
            inputIntervals = readIntervals(bedFile);
        } else {
            inputIntervals = new ArrayList<Interval>();
        }
        for (String r : region) {
            inputIntervals.add(parseRegion(r));
        }

        if (inputIntervals.size() <= 1) return inputIntervals;

        Collections.sort(inputIntervals);

        // Merge overlapping intervals
        ArrayList<Interval> outIntervals = new ArrayList<Interval>(inputIntervals.size());
        Interval activeInterval = outIntervals.get(0);
        for (Interval next : outIntervals) {
            if (activeInterval.abuts(next) || activeInterval.intersects(next)) {
                activeInterval = new Interval(activeInterval.getContig(), activeInterval.getStart(), next.getEnd());
                System.err.printf("WARNING: intervals overlap, coalescing:\n\t%s\n\t%s\n", activeInterval, next);
            } else {
                outIntervals.add(activeInterval);
                activeInterval = next;
            }
        }
        outIntervals.add(activeInterval);

        return outIntervals;
    }

    public void analyze() throws IOException {
        ArrayList<Interval> intervals = setupIntervals();

        for (String bamFile : bamFiles) {
            PrintStream fasta = null;
            PrintStream bedGraph = null;

            String prefix = bamFile;
            if (prefix.endsWith(".bam") || prefix.endsWith(".BAM")) {
                prefix = prefix.substring(0, prefix.length()-4);
            }
            prefix += suffix;

            if (makeCovFasta) {
                if (compressOutput) {
                    fasta = new PrintStream(new BlockCompressedOutputStream(prefix + ".covfasta.gz"));
                } else {
                    fasta = new PrintStream(new BufferedOutputStream(new FileOutputStream(prefix + ".covfasta")));
                }
            }
            if (makeCovBedGraph) {
                if (compressOutput) {
                    bedGraph = new PrintStream(new BlockCompressedOutputStream(prefix + ".bedgraph.gz"));
                } else {
                    bedGraph = new PrintStream(new BufferedOutputStream(new FileOutputStream(prefix + ".bedgraph")));
                }
            }

            PrintStream report = new PrintStream(new File(prefix + ".report"));

            DepthWorker worker = new DepthWorker(minimumMapScore, minimumBaseQuality, keepDupes, bamFile,
                    intervals, fasta, bedGraph, report);
            worker.run();
        }
    }

    public class DepthWorker implements Runnable {
        private final int minimumMapScore;
        private final int minimumBaseQuality;
        private final boolean keepDupes;
        private final String bamFileName;
        private final List<Interval> intervalsOfInterest;
        private final PrintStream covTabOut;
        private final PrintStream covFastaOut;
        private final PrintStream reportOut;

        private final SamReader reader;
        private final long[] coverageHistogram = new long[COVERAGE_HISTOGRAM_MAX+1];

        // maps position -> coverage
        HashMap<Integer, Integer> coverages = new HashMap<Integer, Integer>();


        public DepthWorker(int minimumMapScore, int minimumBaseQuality, boolean keepDupes,
                           String bamFileName, List<Interval> intervalsOfInterest,
                           PrintStream covFastaOut, PrintStream covTabOut, PrintStream reportOut) throws IOException {
            this.minimumMapScore = minimumMapScore;
            this.minimumBaseQuality = minimumBaseQuality;
            this.keepDupes = keepDupes;
            this.bamFileName = bamFileName;

            this.covFastaOut = covFastaOut;
            this.covTabOut = covTabOut;
            this.reportOut = reportOut;

            reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFileName));
            if (intervalsOfInterest.isEmpty()) {
                intervalsOfInterest = wholeGenomeIntervalsForBamFile(reader);
            }
            this.intervalsOfInterest = intervalsOfInterest; // to be really safe, could copy...
        }

        private List<Interval> wholeGenomeIntervalsForBamFile(SamReader reader) {
            ArrayList<Interval> intervals = new ArrayList<Interval>();
            for (SAMSequenceRecord seq : reader.getFileHeader().getSequenceDictionary().getSequences()) {
                intervals.add(new Interval(seq.getSequenceName(), 1, seq.getSequenceLength()));
            }
            return intervals;
        }

        // if this were not java, would closure these currCoverage and currCoverage start into here
        private int currCoverage = -1;
        private Integer currCoverageStart = null;
        private int flushPrior(int alreadyFlushed, int barrier, Interval region) {
            for (int i = alreadyFlushed + 1; i < barrier; i++) {
                Integer coverage = coverages.remove(i);
                if (coverage == null) coverage = 0;
                if (coverage > COVERAGE_HISTOGRAM_MAX) {
                    coverageHistogram[COVERAGE_HISTOGRAM_MAX]++;
                } else {
                    coverageHistogram[coverage]++;
                }
                if (covFastaOut != null) {
                    covFastaOut.print((i-region.getStart()) % 100 == 0 ? "\n" : " ");
                    covFastaOut.print(coverage);
                }
                if (coverage != currCoverage) {
                    if (currCoverageStart != null) {
                        if (covTabOut != null)
                            covTabOut.printf("%s\t%d\t%d\t%d\n", region.getContig(), currCoverageStart-1, i - 1, currCoverage);
                    }
                    currCoverage = coverage;
                    currCoverageStart = i;
                }
            }
            if (barrier >= region.getEnd() + 1) { // finishing region
                if (covTabOut != null)
                    covTabOut.printf("%s\t%d\t%d\t%d\n", region.getContig(), currCoverageStart-1, barrier-1, currCoverage);
                if (covFastaOut != null)
                    covFastaOut.print("\n");
                alreadyFlushed = -1;
                currCoverage = -1;
                currCoverageStart = null;
            } else {
                alreadyFlushed = barrier - 1;
            }
            return alreadyFlushed;
        }

        @Override
        public void run() {

            long totalReadsPaired = 0;
            long totalPairedReadsWithMappedMates = 0;
            long totalAlignedBases = 0;
            long duplicateReads = 0;

            for (Interval interval : intervalsOfInterest) {
                SAMRecordIterator query = reader.query(interval.getContig(), interval.getStart(), interval.getEnd(), false);
                if (covFastaOut != null)
                    covFastaOut.printf(">%s:%d-%d", interval.getContig(), interval.getStart(), interval.getEnd());
                int flushed = interval.getStart()-1;
                while (query.hasNext()) {
                    SAMRecord rec = query.next();
                    if (rec.getMappingQuality() < minimumMapScore || rec.getReadFailsVendorQualityCheckFlag() || rec.getNotPrimaryAlignmentFlag() || rec.getReadUnmappedFlag()) {
                        continue;
                    }
                    if (rec.getReadPairedFlag()) {
                        totalReadsPaired++;
                        if (!rec.getMateUnmappedFlag()) totalPairedReadsWithMappedMates++;
                    }
                    if (rec.getDuplicateReadFlag()) {
                        duplicateReads++;
                        if (!keepDupes) continue;
                    }

                    byte[] baseQualities = rec.getBaseQualities();
                    for (int i = 0; i < rec.getReadLength(); i++) {
                        int refPos = rec.getReferencePositionAtReadPosition(i+1);  // htsjdk expects ONE-based read offset
                        if (refPos > 0) { // 0 is used for 'no corresponding position', e.g. insertion
                            totalAlignedBases++;
                            if (baseQualities[i] >= minimumBaseQuality) coverages.put(refPos, coverages.getOrDefault(refPos, 0) + 1);
                        }
                    }
                    flushed = flushPrior(flushed, rec.getAlignmentStart(), interval);
                }
                flushed = flushPrior(flushed, interval.getEnd()+1, interval);
                query.close();
            }
            try {
                reader.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            if (covFastaOut != null) {
                covFastaOut.close();
            }
            if (covTabOut != null) {
                covTabOut.close();
            }

            //----------------------------   Write report
            reportOut.printf("Total Reads Paired:\t%d\n", totalReadsPaired);
            reportOut.printf("Total Paired Reads With Mapped Mates:\t%d\n", totalPairedReadsWithMappedMates);
            reportOut.printf("Duplicate Reads:\t%d\n", duplicateReads);
            reportOut.printf("Total aligned bases:\t%d\n", totalAlignedBases);
            long hist_at_least[] = new long[coverageHistogram.length];
            long cumulative = 0;
            for (int i = COVERAGE_HISTOGRAM_MAX; i >= 0 ; i--) {
                cumulative += coverageHistogram[i];
                hist_at_least[i] = cumulative;
            }
            cumulative = 0;
            reportOut.println("#-----------------------------------------------------------------------------------------#");
            reportOut.println("Coverage\tCount\tCumulative Below\tCumulativeAbove");
            for (int i = 0; i <= COVERAGE_HISTOGRAM_MAX; i++) {
                cumulative += coverageHistogram[i];
                reportOut.printf("%d\t%d\t%d\t%d\n", i, coverageHistogram[i], cumulative, hist_at_least[i]);
            }
            reportOut.println("#-----------------------------------------------------------------------------------------#");
            reportOut.close();
        }
    }

    private static ArrayList<Interval> readIntervals(String bedFileName) throws IOException {
        ArrayList<Interval> intervals = new ArrayList<Interval>();

        // Using StartOffset.ONE converts to one-based on input
        FeatureReader<BEDFeature> bedReader = AbstractFeatureReader.getFeatureReader(bedFileName, new BEDCodec(BEDCodec.StartOffset.ONE), false);
        BEDFeature f;
        CloseableTribbleIterator<BEDFeature> bedIt = bedReader.iterator();
        while (bedIt.hasNext()) {
            f = bedIt.next();
            intervals.add(new Interval(f.getContig(), f.getStart(), f.getEnd()));
        }
        bedReader.close();
        return intervals;
    }

    public static String joinInts(String separator, int... ints) {
        if (ints.length == 0) return "";
        if (ints.length == 1) return Integer.toString(ints[0]);
        StringBuilder sb = new StringBuilder().append(ints[0]);
        for (int i = 1; i < ints.length; i++) {
            sb.append(separator).append(ints[i]);
        }
        return sb.toString();
    }

    public static PrintStream openOutput(String filename) throws IOException {
        if (filename.endsWith(".gz")) {
            return new PrintStream(new GZIPOutputStream(new FileOutputStream(filename)));
        } else {
            return new PrintStream(new BufferedOutputStream(new FileOutputStream(filename)));
        }
    }

    private static Interval parseRegion(String region) {
        int colon = region.lastIndexOf(':');
        String contig = region.substring(0, colon);
        int dash = region.lastIndexOf('-');
        int start = Integer.parseInt(region.substring(colon+1, dash));
        int end = Integer.parseInt(region.substring(dash+1));
        return new Interval(contig, start, end);
    }

    public static void main(String [] args) throws IOException {
        CallableDepth cd = new CallableDepth();
        JCommander jc = new JCommander(cd);
        jc.parse(args);
        if (cd.help) {
            jc.usage();
            System.exit(0);
        }
        if (cd.bamFiles.size() != 1) {
            jc.usage();
            throw new IllegalArgumentException("Must provide exactly one bam file!");
        }

        if (!cd.suffix.isEmpty() && !cd.suffix.startsWith(".")) {
            cd.suffix = "." + cd.suffix;
        }

        cd.analyze();
    }

}
