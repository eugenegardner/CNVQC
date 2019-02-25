package merger.cluster;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import utilities.CNV;
import utilities.CNVInterval;

public interface CNVClusterer {
		
	public Map<Integer, CNVInterval> getMergedCNVs(List<CNV> currentCNVList) throws IOException;
	public ClusterType getClusterType();
	
	public enum ClusterType {
		DBSCAN,RECIPROCAL;
	}
	
}
