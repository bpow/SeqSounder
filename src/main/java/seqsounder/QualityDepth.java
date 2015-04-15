package seqsounder;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import seqsounder.depthresponder.*;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import java.io.*;
import java.lang.reflect.Field;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class QualityDepth {
    private static final int COVERAGE_HISTOGRAM_MAX = 1000;
    private int minimumMapScore = 20;
    private int minimumBaseQuality = 20;
    private boolean keepDupes = false;
    private boolean compressOutput = true;
    private String bedFile = null;
    private boolean makeCovFasta = false;
    private boolean makeCovBedGraph = true;
    private String suffix = "";
    private List<String> bamFiles = new ArrayList<String>();
    private int threadCount = 1;

    public QualityDepth set(String parameter, Object value) {
        Class<? extends QualityDepth> clazz = this.getClass();
        try {
            Field f = clazz.getDeclaredField(parameter);
            f.set(this, value);
        } catch (NoSuchFieldException e) {
            e.printStackTrace();
            return null;
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        }
        return this;
    }

    public void analyze() throws IOException {
        ArrayList<Interval> intervals = bedFile == null ? new ArrayList<Interval>() : readIntervals(bedFile);
        final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

        for (String bamFile : bamFiles) {
            PrintStream fasta = null;
            PrintStream bedGraph = null;

            String prefix = bamFile;
            if (prefix.endsWith(".bam") || prefix.endsWith(".BAM")) {
                prefix = prefix.substring(0, prefix.length()-4);
            }
            prefix += suffix;

            AggregatingResponder aggregator = new AggregatingResponder();

            if (makeCovFasta) {
                if (compressOutput) {
                    fasta = new PrintStream(new BlockCompressedOutputStream(prefix + ".covfasta.gz"));
                } else {
                    fasta = new PrintStream(new BufferedOutputStream(new FileOutputStream(prefix + ".covfasta")));
                }
                aggregator.addClients(new CovFastaResponder(fasta));
            }
            if (makeCovBedGraph) {
                if (compressOutput) {
                    bedGraph = new PrintStream(new BlockCompressedOutputStream(prefix + ".bedgraph.gz"));
                } else {
                    bedGraph = new PrintStream(new BufferedOutputStream(new FileOutputStream(prefix + ".bedgraph")));
                }
                aggregator.addClients(new BedGraphResponder(bedGraph));
            }

            DepthWorker worker = new DepthWorker(minimumMapScore, minimumBaseQuality, keepDupes, COVERAGE_HISTOGRAM_MAX, bamFile,
                    intervals, aggregator);
            executor.execute(worker);
        }
        try {
            executor.awaitTermination(1024L, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        executor.shutdown();
    }

    public static ArrayList<Interval> readIntervals(String bedFileName) throws IOException {
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
        if (intervals.size() <= 1) return intervals;

        Collections.sort(intervals);

        // Merge overlapping intervals
        ArrayList<Interval> checkedIntervals = new ArrayList<Interval>(intervals.size());
        Interval activeInterval = intervals.remove(0);
        for (Interval next : intervals) {
            if (activeInterval.abuts(next) || activeInterval.intersects(next)) {
                activeInterval = new Interval(activeInterval.getContig(), activeInterval.getStart(), next.getEnd());
                System.err.printf("WARNING: intervals overlap, coalescing:\n\t%s\n\t%s\n", activeInterval, next);
            } else {
                checkedIntervals.add(activeInterval);
                activeInterval = next;
            }
        }
        checkedIntervals.add(activeInterval);

        return checkedIntervals;

    }

    public QualityDepth applyConfig(InputStream config) throws ScriptException {
        ScriptEngine engine = new ScriptEngineManager().getEngineByName("JavaScript");
        Scanner s = new Scanner(config);
        String jsString = s.useDelimiter("\\A").next(); // slurp whole file
        s.close();
        engine.put("qd", this);
        engine.eval(jsString);
        return this;
    }

    public static void main(String [] args) throws IOException, ScriptException {
        QualityDepth qd = new QualityDepth();
        qd.applyConfig(new FileInputStream(args[0]));
        qd.analyze();
    }

}
