package seqsounder;

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
import seqsounder.depthresponder.*;

import java.io.*;
import java.util.*;

public class QualityDepth {
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
        Interval activeInterval = inputIntervals.remove(0);
        for (Interval next : inputIntervals) {
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

    private class CoverageMap extends HashMap<Integer, SiteDepth> {
        public SiteDepth getDefault(Integer key) {
            SiteDepth mi = super.get(key);
            if (null == mi) {
                mi = new SiteDepth(key, 0);
                super.put(key, mi);
            }
            return mi;
        }
    }

    public class DepthWorker implements Runnable {
        private final int minimumMapScore;
        private final int minimumBaseQuality;
        private final boolean keepDupes;
        private final List<Interval> intervalsOfInterest;
        private final PrintStream covTabOut;
        private final PrintStream covFastaOut;
        private final PrintStream reportOut;
        private final ArrayList<DepthResponder> responders = new ArrayList<DepthResponder>();

        private final SamReader reader;
        private final HistogramResponder histogramCreator = new HistogramResponder(COVERAGE_HISTOGRAM_MAX);

        // maps position -> coverage
        CoverageMap coverages = new CoverageMap();

        public DepthWorker(int minimumMapScore, int minimumBaseQuality, boolean keepDupes,
                           String bamFileName, List<Interval> intervalsOfInterest,
                           PrintStream covFastaOut, PrintStream covTabOut, PrintStream reportOut) throws IOException {
            this.minimumMapScore = minimumMapScore;
            this.minimumBaseQuality = minimumBaseQuality;
            this.keepDupes = keepDupes;

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

        private int flushPrior(int alreadyFlushed, int barrier, Interval region) {
            if (alreadyFlushed + 1 >= barrier) return alreadyFlushed;
            for (int i = alreadyFlushed + 1; i < barrier; i++) {
                SiteDepth mutCoverage = coverages.remove(i);
                if (mutCoverage == null) {
                    mutCoverage = new SiteDepth(i, 0);
                }
                for (DepthResponder responder : responders) {
                    responder.markDepth(mutCoverage);
                }
            }
            alreadyFlushed = barrier - 1;
            if (barrier >= region.getEnd() + 1) { // finishing region
                alreadyFlushed = -1;
            } else {
            }
            return alreadyFlushed;
        }

        @Override
        public void run() {

            long totalReadsPaired = 0;
            long totalPairedReadsWithMappedMates = 0;
            long totalAlignedBases = 0;
            long duplicateReads = 0;

            AggregatingResponder aggregator = new AggregatingResponder(histogramCreator);
            if (covTabOut != null) {
                aggregator.addClients(new BedGraphResponder(covTabOut));
            }
            if (covFastaOut != null) {
                aggregator.addClients(new CovFastaResponder(covFastaOut));
            }
            responders.add(aggregator);

            for (Interval interval : intervalsOfInterest) {
                SAMRecordIterator query = reader.query(interval.getContig(), interval.getStart(), interval.getEnd(), false);
                for (DepthResponder responder : responders) {
                    responder.startRegion(interval);
                }
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
                        if (refPos >= interval.getStart() && refPos <= interval.getEnd()) { // 0 is used for 'no corresponding position', e.g. insertion
                            totalAlignedBases++;
                            if (baseQualities[i] >= minimumBaseQuality) coverages.getDefault(refPos).increment();
                        }
                    }
                    flushed = flushPrior(flushed, rec.getAlignmentStart(), interval);
                }
                flushed = flushPrior(flushed, interval.getEnd()+1, interval);
                for (DepthResponder responder : responders) {
                    responder.finishRegion(interval);
                }
                query.close();
            }
            for (DepthResponder responder : responders) {
                responder.finishAll();
            }
            try {
                reader.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            //----------------------------   Write report
            reportOut.printf("Total Reads Paired:\t%d\n", totalReadsPaired);
            reportOut.printf("Total Paired Reads With Mapped Mates:\t%d\n", totalPairedReadsWithMappedMates);
            reportOut.printf("Duplicate Reads:\t%d\n", duplicateReads);
            reportOut.printf("Total aligned bases:\t%d\n", totalAlignedBases);
            long [] coverageHistogram = histogramCreator.getHistogram();
            long hist_at_least[] = new long[coverageHistogram.length];
            long cumulative = 0;
            for (int i = coverageHistogram.length-1; i >= 0 ; i--) {
                cumulative += coverageHistogram[i];
                hist_at_least[i] = cumulative;
            }
            cumulative = 0;
            reportOut.println("#-----------------------------------------------------------------------------------------#");
            reportOut.println("Coverage\tCount\tCumulative Below\tCumulative Above");
            for (int i = 0; i < coverageHistogram.length; i++) {
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

    private static Interval parseRegion(String region) {
        int colon = region.lastIndexOf(':');
        String contig = region.substring(0, colon);
        int dash = region.lastIndexOf('-');
        int start = Integer.parseInt(region.substring(colon+1, dash));
        int end = Integer.parseInt(region.substring(dash+1));
        return new Interval(contig, start, end);
    }

    public static void main(String [] args) throws IOException {
        QualityDepth qd = new QualityDepth();
        JCommander jc = new JCommander(qd);
        jc.parse(args);
        if (qd.help) {
            jc.usage();
            System.exit(0);
        }
        if (qd.bamFiles.size() != 1) {
            jc.usage();
            throw new IllegalArgumentException("Must provide exactly one bam file!");
        }

        if (!qd.suffix.isEmpty() && !qd.suffix.startsWith(".")) {
            qd.suffix = "." + qd.suffix;
        }

        qd.analyze();
    }

}
