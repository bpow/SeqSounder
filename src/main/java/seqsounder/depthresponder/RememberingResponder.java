package seqsounder.depthresponder;

import seqsounder.SiteDepth;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class RememberingResponder extends DepthResponder {
    private final ArrayList<SiteDepth> depths = new ArrayList<SiteDepth>();

    @Override
    public void markDepth(SiteDepth siteDepth) {
        depths.add(siteDepth);
    }

    public List<SiteDepth> getDepths() {
        return Collections.unmodifiableList(depths);
    }
}
