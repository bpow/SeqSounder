package seqsounder.depthresponder;

import htsjdk.samtools.util.Interval;
import seqsounder.SiteDepth;

public abstract class DepthResponder {
    public void startRegion(Interval region) {}
    public abstract void markDepth(SiteDepth siteDepth);
    public void finishRegion(Interval region) {}
    public void finishAll() {}
}
