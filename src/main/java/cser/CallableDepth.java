package cser;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPOutputStream;

public class CallableDepth {
    private final int minimumMapScore;
    private final int minimumBaseQuality;
    private final boolean keepDupes;
    private final String bamFileName;
    private final List<Interval> queryIntervals;
    private final PrintStream covFastaOut;
    private final PrintStream covTabOut;

    public CallableDepth(int minimumMapScore, int minimumBaseQuality, boolean keepDupes,
                         String bamFileName, List<Interval> queryIntervals,
                         String outCovFasta, String outCovTab) throws FileNotFoundException {
        this.minimumMapScore = minimumMapScore;
        this.minimumBaseQuality = minimumBaseQuality;
        this.keepDupes = keepDupes;
        this.bamFileName = bamFileName;
        this.queryIntervals = queryIntervals;

        // setup outputs
        if (outCovTab == null) {
            covTabOut = null;
            if (outCovFasta == null) {
                covFastaOut = System.out; // have to have some output...
            } else {
                covFastaOut = new PrintStream(outCovFasta);
            }
        } else {
            covTabOut = new PrintStream(outCovTab);
            covFastaOut = outCovFasta == null ? null : new PrintStream(outCovFasta);
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

    public void analyze() {
        SamReaderFactory srf = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT);
        SamReader reader = srf.open(new File(bamFileName));

        long totalReadsPaired = 0;
        long totalPairedReadWithMappedMates = 0;
        long totalAlignedBases = 0;
        long duplicateReads = 0;

        for (Interval interval : queryIntervals) {
            int[] coverages = new int[interval.length()];
            SAMRecordIterator query = reader.queryOverlapping(interval.getContig(), interval.getStart(), interval.getEnd());
            // cannot just use enhanced for loop (or groovy 'each') because htsjdk needs its iterators to be `close`d, so we need a reference to the iterator!
            while (query.hasNext()) {
                SAMRecord rec = query.next();

                if (rec.getMappingQuality() < minimumMapScore || rec.getReadFailsVendorQualityCheckFlag() || rec.getNotPrimaryAlignmentFlag() || rec.getReadUnmappedFlag()) {
                    continue;
                }
                if (rec.getReadPairedFlag()) {
                    totalReadsPaired++;
                    if (!rec.getMateUnmappedFlag()) totalPairedReadWithMappedMates++;
                }
                if (rec.getDuplicateReadFlag()) {
                    duplicateReads++;
                    if (!keepDupes) continue;
                }
                byte[] baseQualities = rec.getBaseQualities();
                for (int i = 0; i < rec.getReadLength(); i++) {
                    int refPos = rec.getReferencePositionAtReadPosition(i + 1); // htsjdk expects ONE-based read offset
                    if (refPos > 0) { // 0 is used for 'no corresponding position', e.g. insertion
                        int intervalPos = refPos - interval.getStart(); // ZERO-based intervalPos
                        if ((intervalPos >= 0) && (intervalPos < interval.length())) {
                            totalAlignedBases++;
                            if (baseQualities[i] >= minimumBaseQuality) coverages[intervalPos]++;
                        }
                    }
                }
            }
            query.close();
            flushOutput(interval, coverages);
        }
        if (covFastaOut != null) {
            covFastaOut.close();
        }
        if (covTabOut != null) {
            covTabOut.close();
        }

    }

    private void flushOutput(Interval interval, int[] coverages) {
        assert coverages.length == (interval.getEnd()-interval.getStart()+1);
        if (covFastaOut != null) {
            covFastaOut.printf(">%s:%d-%d\n", interval.getContig(), interval.getStart(), interval.getEnd());
            covFastaOut.println(joinInts(" ", coverages));
        }
        if (covTabOut != null) {
            int start = interval.getStart();
            String contig = interval.getContig();
            for (int i = 0; i < coverages.length; i++) {
                covTabOut.printf("%s\t%d\t%d", interval.getContig(), i+start, coverages[i]);
            }
        }
    }

    public static String joinInts(String separator, int... ints) {
        if (ints.length == 0) return "";
        if (ints.length == 0) return Integer.toString(ints[0]);
        StringBuilder sb = new StringBuilder().append(ints[0]);
        for (int i = 1; i < ints.length; i++) {
            sb.append(separator).append(ints[i]);
        }
        return sb.toString();
    }

    public static PrintWriter openOutput(String filename) throws IOException {
        if (filename.endsWith(".gz")) {
            return new PrintWriter(new GZIPOutputStream(new FileOutputStream(filename)));
        } else {
            return new PrintWriter(filename);
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

    static class ArgParser {
        @Parameter(names = {"-h", "--help"}, help = true)
        private boolean help;
        @Parameter(names="-mapQ", description="Minimum mapping quality")
        public int minimumMapScore = 20;
        @Parameter(names="-baseQ", description="Minimum base quality")
        public int minimumBaseQuality = 20;
        @Parameter(names="-keepDupes", description="include duplicates in analysis")
        public boolean keepDupes;
        @Parameter(names={"-bed", "-l"}, description="BED file of intervals to analyze (ZERO-based coordinates)")
        public String bedFile = null;
        @Parameter(names={"-r", "-region"}, description="Regions to analyze (ONE-based), like \"1:32221-42212\", can be specified multiple times")
        public List<String> region;
        @Parameter(names={"-o", "-covFasta"}, description="Ouput as a \"coverage.fasta\" file (will output to stdout if no output specified)")
        public String covFastaOut = null;
        @Parameter(names={"-t", "-covTSV"}, description="Ouput as a tab-delimited text file")
        public String covTabOut = null;
        @Parameter(description = "bamFile")
        public List<String> bamFiles = new ArrayList<String>();
    }

    public static void main(String [] args) throws IOException {
        ArgParser ap = new ArgParser();
        JCommander jc = new JCommander(ap);
        jc.parse(args);
        if (ap.help) {
            jc.usage();
            System.exit(0);
        }
        if (ap.bamFiles.size() != 1) {
            jc.usage();
            throw new IllegalArgumentException("Must provide exactly one bam file!");
        }

        ArrayList<Interval> queryIntervals;
        if (ap.region.isEmpty() && ap.bedFile == null) {
            jc.usage();
            throw new IllegalArgumentException("Must provide one of '-r' or '-l'");
        } else {
            if (ap.bedFile != null) {
                queryIntervals = CallableDepth.readIntervals(ap.bedFile);
            } else {
                queryIntervals = new ArrayList<Interval>();
            }
            for (String region : ap.region) {
                queryIntervals.add(parseRegion(region));
            }
            Collections.sort(queryIntervals);
        }

        CallableDepth main = new CallableDepth(ap.minimumMapScore, ap.minimumBaseQuality, ap.keepDupes,
                ap.bamFiles.get(0), queryIntervals, ap.covFastaOut, ap.covTabOut);
        main.analyze();
    }

}
