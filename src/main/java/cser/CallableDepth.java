package cser;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeEOFException;
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
    @Parameter(names = {"-h", "--help"}, help = true)
    public boolean help;
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
    public String covFasta = null;
    @Parameter(names={"-t", "-covTab"}, description="Ouput as a tab-delimited text file")
    public String covTab = null;
    @Parameter(description = "bamFile")
    public List<String> bamFiles = new ArrayList<String>();

    public void analyze() throws IOException {
        ArrayList<Interval> queryIntervals;

        if (bedFile != null) {
            queryIntervals = readIntervals(bedFile);
        } else {
            queryIntervals = new ArrayList<Interval>();
        }
        for (String r : region) {
            queryIntervals.add(parseRegion(r));
        }
        Collections.sort(queryIntervals);


        DepthWorker worker = new DepthWorker(minimumMapScore, minimumBaseQuality, keepDupes, bamFiles.get(0),
                queryIntervals, covFasta, covTab);
        worker.run();
    }

    public class DepthWorker implements Runnable {
        private final int minimumMapScore;
        private final int minimumBaseQuality;
        private final boolean keepDupes;
        private final String bamFileName;
        private final List<Interval> queryIntervals;
        private final PrintStream covTabOut;
        private final PrintStream covFastaOut;

        private final SamReader reader;

        // maps position -> coverage
        HashMap<Integer, Integer> coverages = new HashMap<Integer, Integer>();
        String currContig = null;
        int flushed = 0;


        public DepthWorker(int minimumMapScore, int minimumBaseQuality, boolean keepDupes,
                           String bamFileName, List<Interval> queryIntervals,
                           String outCovFasta, String outCovTab) throws FileNotFoundException {
            this.minimumMapScore = minimumMapScore;
            this.minimumBaseQuality = minimumBaseQuality;
            this.keepDupes = keepDupes;
            this.bamFileName = bamFileName;
            this.queryIntervals = queryIntervals; // to be really safe, could copy...

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
            reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFileName));
        }

        private void flushPrior(int barrier) {
            for (int i = flushed + 1; i < barrier; i++) {
                Integer coverage = coverages.remove(i);
                if (coverage == null) coverage = 0;
                covTabOut.printf("%s\t%d\t%d\n", currContig, i, coverage);
                flushed = i;
            }
        }

        @Override
        public void run() {

            long totalReadsPaired = 0;
            long totalPairedReadWithMappedMates = 0;
            long totalAlignedBases = 0;
            long duplicateReads = 0;

            for (SAMRecord rec : reader) {
                if (rec.getContig() != currContig && currContig != null) {
                    flushPrior(reader.getFileHeader().getSequence(currContig).getSequenceLength()+1);
                }
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
                    int refPos = rec.getReferencePositionAtReadPosition(i+1);  // htsjdk expects ONE-based read offset
                    if (refPos > 0) { // 0 is used for 'no corresponding position', e.g. insertion
                        totalAlignedBases++;
                        if (baseQualities[i] >= minimumBaseQuality) coverages.put(refPos, coverages.getOrDefault(refPos, 0) + 1);
                    }
                }
                currContig = rec.getContig();
                flushPrior(rec.getAlignmentStart());

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

        ArrayList<Interval> queryIntervals;
        if (cd.region.isEmpty() && cd.bedFile == null) {
            jc.usage();
            throw new IllegalArgumentException("Must provide one of '-r' or '-l'");
        }

        cd.analyze();
    }

}
