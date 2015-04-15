package seqsounder.depthresponder;

import htsjdk.samtools.util.Interval;
import seqsounder.SiteDepth;

import java.util.ArrayList;
import java.util.Arrays;

public class AggregatingResponder extends DepthResponder {
    private SiteDepth workingSiteDepth = null;
    private final ArrayList<DepthResponder> clients = new ArrayList<DepthResponder>();

    @Override
    public void startRegion(Interval region) {
        for (DepthResponder client : clients) {
            client.startRegion(region);
        }
        workingSiteDepth = null;
    }

    @Override
    public void markDepth(SiteDepth siteDepth) {
        if (workingSiteDepth == null) {
            workingSiteDepth = siteDepth;
        } else {
            if (workingSiteDepth.end < siteDepth.start) { // discontinuity
                if (workingSiteDepth.depth == 0) {
                    workingSiteDepth = new SiteDepth(workingSiteDepth.start, siteDepth.start, workingSiteDepth.depth);
                } else {
                    emitDepth(workingSiteDepth);
                    workingSiteDepth = new SiteDepth(workingSiteDepth.end, siteDepth.start, 0);
                }
            }
            if (workingSiteDepth.depth != siteDepth.depth) {
                emitDepth(workingSiteDepth);
                workingSiteDepth = siteDepth;
            } else {
                workingSiteDepth = new SiteDepth(workingSiteDepth.start, siteDepth.end, siteDepth.depth);
            }
        }
    }

    @Override
    public void finishRegion(Interval region) {
        if (workingSiteDepth != null) {
            emitDepth(workingSiteDepth);
            workingSiteDepth = null;
        }
        for (DepthResponder client : clients) {
            client.finishRegion(region);
        }
    }

    @Override
    public void finishAll() {
        for (DepthResponder client : clients) {
            client.finishAll();
        }
    }

    private final void emitDepth(SiteDepth siteDepth) {
        for (DepthResponder client : clients) {
            client.markDepth(siteDepth);
        }
    }

    public AggregatingResponder(DepthResponder... clients) {
        this.clients.addAll(Arrays.asList(clients));
    }

    public AggregatingResponder addClients(DepthResponder... clients) {
        this.clients.addAll(Arrays.asList(clients));
        return this;
    }
}
