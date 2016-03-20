package seqsounder.irods;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.irods.jargon.core.connection.IRODSAccount;
import org.irods.jargon.core.exception.JargonException;
import org.irods.jargon.core.pub.IRODSFileSystem;
import org.irods.jargon.core.pub.domain.DataObject;
import org.irods.jargon.core.pub.io.IRODSFile;
import org.irods.jargon.core.pub.io.IRODSFileFactory;
import org.irods.jargon.core.query.AVUQueryElement;
import org.irods.jargon.core.query.AVUQueryOperatorEnum;
import org.irods.jargon.core.query.JargonQueryException;
import org.irods.jargon.core.utils.IRODSUriUtils;

import java.net.URI;
import java.util.*;

public class IrodsQueryUtils {
    private final IRODSFileSystem fs;
    private final IRODSAccount account;

    public List<DataObject> irodsAVQuery(Map<String, String> keyValues) throws JargonQueryException, JargonException {
        ArrayList<AVUQueryElement> queryElements = new ArrayList<AVUQueryElement>(2*keyValues.size());
        for (Map.Entry<String, String> entry: keyValues.entrySet()) {
            queryElements.add(AVUQueryElement.instanceForValueQuery(
                    AVUQueryElement.AVUQueryPart.ATTRIBUTE,
                    AVUQueryOperatorEnum.EQUAL, entry.getKey()));
            queryElements.add(AVUQueryElement.instanceForValueQuery(
                    AVUQueryElement.AVUQueryPart.VALUE,
                    AVUQueryOperatorEnum.EQUAL, entry.getValue()));
        }
        List<DataObject> results = fs.getIRODSAccessObjectFactory()
                .getDataObjectAO(account).findDomainByMetadataQuery(queryElements);
        results.sort(new Comparator<DataObject>() {
            @Override
            public int compare(DataObject o1, DataObject o2) {
                return o2.getCreatedAt().compareTo(o1.getCreatedAt());
            }
        });
        return results;
    }

    private IrodsSeekableStream streamForDataObject(DataObject dataObj) throws JargonException {
        IRODSFileFactory irff = fs.getIRODSFileFactory(account);
        IRODSFile irf = irff.instanceIRODSFile(dataObj.getAbsolutePath());
        irf.openReadOnly();
        // in jargon 4.0, will be:
        // irf.open(DataObjInp.OpenFlags.READ);
        return new IrodsSeekableStream(irff.instanceIRODSRandomAccessFile(irf));
    }

    public SamReader readerForBam(String bamFileName) {
        final String participantID = bamFileName; // FIXME -- handle '.bam' suffix?
        try {
            List<DataObject> bams = irodsAVQuery(pairsToMap(
                    "ParticipantID", participantID, "FileType", "RecalBam"));
            List<DataObject> bais = irodsAVQuery(pairsToMap(
                    "ParticipantID", participantID, "FileType", "RecalBamBai"));
            IrodsSeekableStream bamStream = streamForDataObject(bams.get(0));
            IrodsSeekableStream baiStream = streamForDataObject(bais.get(0));
            SamInputResource sir = SamInputResource.of(bamStream).index(baiStream);
            SamReader reader = SamReaderFactory.makeDefault().open(sir);
            return reader;
        } catch (JargonQueryException e) {
            e.printStackTrace();
        } catch (JargonException e) {
            e.printStackTrace();
        }
        return null;
    }

    // sooo much easier in groovy, or python, or javascript...
    public static final Map<String, String> pairsToMap(String... pairs) {
        if (pairs.length % 2 != 0) {
            throw new IllegalArgumentException("expected an even number of arguments");
        }
        HashMap<String, String> m = new LinkedHashMap<String, String>(4 * pairs.length / 3);
        for (int i = 0; i < pairs.length; i += 2) {
            m.put(pairs[i], pairs[i+1]);
        }
        return m;
    }

    public IrodsQueryUtils(String irodsUri) throws JargonException {
        fs = new IRODSFileSystem();
        IRODSAccount tempAccount = IRODSUriUtils.getIRODSAccountFromURI(URI.create(irodsUri));
        account = fs.getIRODSAccessObjectFactory().authenticateIRODSAccount(tempAccount).getAuthenticatedIRODSAccount();
    }
}
