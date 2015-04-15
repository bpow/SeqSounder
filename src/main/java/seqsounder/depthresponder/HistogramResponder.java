package seqsounder.depthresponder;

import htsjdk.samtools.util.Interval;
import seqsounder.SiteDepth;

public class HistogramResponder extends DepthResponder {
    private final long [] histogram;

    public HistogramResponder(int maxCoverage) {
        this.histogram = new long[maxCoverage + 1];
    }

    @Override
    public void markDepth(SiteDepth siteDepth) {
        int depth = siteDepth.depth;
        if (depth >= histogram.length) depth = histogram.length-1;
        histogram[depth] += siteDepth.end - siteDepth.start;
    }

    public long [] getHistogram() {
        return histogram;
    }
}
