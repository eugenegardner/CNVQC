package utilities;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import merger.CNVConverter.CopyType;
import sampleannotator.RawCNVReader.LRRandBAFInformation;
import sampleannotator.resources.SampleInformation;

public class CNV {

	private String chr;
	private int start;
	private int end;
	private int copyNumber;
	private int probeCount;
	private double confidence;
	private boolean siteFiltered;
	private SampleInformation sampleInformation;
	private LRRandBAFInformation lrrbaf;
	private int mergeGroup;
	private DescriptiveStatistics WESstats;
	private double intersectingWESConvexCNVs;
	private double intersectingWESXHMMCNVs;
	private double intersectingWESCLAMMSCNVs;
	private double intersectingWESCANOESCNVs;
	private String gsCNV;
	private double distTel;
	private double distCen;
	private String printable;

	/**
	 * @param chr Chromosome on which CNV resides
	 * @param start Start coordinate of CNV
	 * @param end End coordinate of CNV
	 * @param copyNumber Integer copy number value
	 * @param probeCount Total number of overlapping Affy probes
	 * @param confidence PennCNV confidence score (equals probability of highest copy number state minus probability of next most likely state)
	 * @param sampleInformation Collection of information pertaining to the sample in which this CNV was found
	 * @param lrrbaf Collection of LRR and BAF scores surrounding the CNV
	 */
	public CNV (String chr, int start, int end, int copyNumber, int probeCount, double confidence, SampleInformation sampleInformation, LRRandBAFInformation lrrbaf) {
		
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.copyNumber = copyNumber;
		this.probeCount = probeCount;
		this.confidence = confidence;
		this.sampleInformation = sampleInformation;
		this.lrrbaf = lrrbaf;
		
		siteFiltered = checkSiteFilter();
		
		totalIntersectingBaits = 0;
		WESstats = new DescriptiveStatistics();
		
		mergeGroup = -1;
		
		this.printable = null;
		
	}
	public CNV (String chr, int start, int end, int copyNumber, int probeCount, double confidence, SampleInformation sampleInformation, LRRandBAFInformation lrrbaf, String printable) {
		
		this(chr, start, end, copyNumber, probeCount, confidence, sampleInformation, lrrbaf);
		this.printable = printable;
		
	}
	

	public void attachMergeGroup(int mergeGroup) {
		this.mergeGroup = mergeGroup;
	}
	public int getMergeGroup() {
		return mergeGroup;
	}
	
	public int getLength() {
		return end - start;
	}
	public String getChr() {
		return chr;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getCopyNumber() {
		return copyNumber;
	}
	public String getLocationCoordinates() {
		return chr + ":" + start + "-" + end;
	}	
	public int getProbeCount() {
		return probeCount;
	}
	public double getConfidence() {
		return confidence;
	}
	public double getDensity() {
		double density = (double) this.getLength() / (double) probeCount;
		return density;
	}
	public boolean isSiteFiltered() {
		return siteFiltered;
	}
	public SampleInformation getSampleInformation() {
		return sampleInformation;
	}
	public LRRandBAFInformation getLRRBAF() {
		return lrrbaf;
	}
	
	public CopyType getCopyType() {
		if (this.getCopyNumber() < 2) {
			return CopyType.DEL;
		} else {
			return CopyType.DUP;
		}
	}

	private boolean checkSiteFilter() {
		if (this.getDensity() < 20000 && probeCount > 10) {
			return false;
		} else {
			return true;
		}
	}
		
	public void setIntersectingWESConvexCNVs(double intersectingWESConvexCNVs) {
		this.intersectingWESConvexCNVs = intersectingWESConvexCNVs;
	}
	public double getIntersectingWESConvexCNVs() {
		return intersectingWESConvexCNVs;
	}
	public void setIntersectingWESXHMMCNVs(double intersectingWESXHMMCNVs) {
		this.intersectingWESXHMMCNVs = intersectingWESXHMMCNVs;
	}
	public double getIntersectingWESXHMMCNVs() {
		return intersectingWESXHMMCNVs;
	}
	public void setIntersectingWESCLAMMSCNVs(double intersectingWESCLAMMSCNVs) {
		this.intersectingWESCLAMMSCNVs = intersectingWESCLAMMSCNVs;
	}
	public double getIntersectingWESCLAMMSCNVs() {
		return intersectingWESCLAMMSCNVs;
	}
	public void setIntersectingWESCANOESCNVs(double intersectingWESCANOESCNVs) {
		this.intersectingWESCANOESCNVs = intersectingWESCANOESCNVs;
	}
	public double getIntersectingWESCANOESCNVs() {
		return intersectingWESCANOESCNVs;
	}
		
	private int totalIntersectingBaits;
	
	public int getTotalIntersectingBaits() {
		return totalIntersectingBaits;
	}
	public void setTotalIntersectingBaits(int totalIntersectingBaits) {
		this.totalIntersectingBaits = totalIntersectingBaits;
	}
	public void setWESStats(DescriptiveStatistics WESstats) {
		this.WESstats = WESstats;
	}
	public double getWESMeanL2R() {
		return WESstats.getMean();
	}
	public double getWESL2Rsd() {
		return WESstats.getStandardDeviation();
	}
	public long getWESL2N() {
		return WESstats.getN();
	}

	public void setGoldStandardSV(String gsCNV) {
		this.gsCNV = gsCNV;
	}
	public String getGoldStandardSV() {
		return gsCNV;
	}
	public double getDistTel() {
		return distTel;
	}
	public double getDistCen() {
		return distCen;
	}
	public void setDistTel(double d) {
		this.distTel = d;
	}
	public void setDistCen(double d) {
		this.distCen = d;
	}
	
	public String getPrintable() {
		return printable;
	}


}
