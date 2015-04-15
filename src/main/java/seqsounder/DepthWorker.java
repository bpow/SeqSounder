package seqsounder;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import seqsounder.depthresponder.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

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
    private final HistogramResponder histogramCreator;

    // maps position -> coverage
    CoverageMap coverages = new CoverageMap();

    public DepthWorker(int minimumMapScore, int minimumBaseQuality, boolean keepDupes, int histogramMaxDepth,
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
        histogramCreator = new HistogramResponder(histogramMaxDepth);
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
            int flushed = interval.getStart() - 1;
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
                    int refPos = rec.getReferencePositionAtReadPosition(i + 1);  // htsjdk expects ONE-based read offset
                    if (refPos >= interval.getStart() && refPos <= interval.getEnd()) { // 0 is used for 'no corresponding position', e.g. insertion
                        totalAlignedBases++;
                        if (baseQualities[i] >= minimumBaseQuality) coverages.getDefault(refPos).increment();
                    }
                }
                flushed = flushPrior(flushed, rec.getAlignmentStart(), interval);
            }
            flushed = flushPrior(flushed, interval.getEnd() + 1, interval);
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
        long[] coverageHistogram = histogramCreator.getHistogram();
        long hist_at_least[] = new long[coverageHistogram.length];
        long cumulative = 0;
        for (int i = coverageHistogram.length - 1; i >= 0; i--) {
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


}
