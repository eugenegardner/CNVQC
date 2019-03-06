package sampleannotator;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import merger.CNVConverter.CopyType;
import sampleannotator.resources.SampleInformation;
import utilities.CNV;

public class RawCNVReader implements Closeable {

	private Set<String> samples;
	private Set<String> chrs;
	private Map<String, SampleInformation> sampleInformation;
	private BufferedReader cnvReader;
	private TabixReader axiomTabix;
	private TabixReader baitTabix;
	private Map<String,IntervalTree<String>> dgvDels;
	private Map<String,IntervalTree<String>> dgvDups;
	Map<String, CytobandInfo> cytoMap;
	
	public RawCNVReader(File CNVs, Map<String, SampleInformation> sampleInformation) throws IOException {
		
		cnvReader = new BufferedReader(new FileReader(CNVs));
		samples = new HashSet<String>();
		chrs = new HashSet<String>();
		this.sampleInformation = sampleInformation;
				
		// Build GS Del/Dups to parse against our CNV set
		dgvDels = new HashMap<String,IntervalTree<String>>();
		dgvDups = new HashMap<String,IntervalTree<String>>();
		BuildGoldStandard();
		
		File axiomProbes = new File("/gluster/neurocluster/projects/biobank-cnv/BATCHES/CNVQC/referencefiles/axiom_probes.bed.gz");
		axiomTabix = new TabixReader(axiomProbes.getAbsolutePath(),axiomProbes.getAbsolutePath() + ".tbi");
		
		File baits = new File("/gluster/neurocluster/projects/biobank-cnv/BATCHES/CNVQC/referencefiles/bait_regions_unpadded.merge500.bed.gz");
		baitTabix = new TabixReader(baits.getAbsolutePath(),baits.getAbsolutePath() + ".tbi");
		
		cytoMap = buildCytoband(new File("/gluster/neurocluster/projects/biobank-cnv/BATCHES/CNVQC/referencefiles/cytoBand.txt"));
		
	}
	
	public Set<String> getChrs() {
		return chrs;
	}
	public Set<String> getSamples() {
		return samples;
	}
	public void close() throws IOException {
		cnvReader.close();
		baitTabix.close();
		axiomTabix.close();

	}
	
	private void BuildGoldStandard() throws IOException {
		File CNVGoldStandard = new File("/gluster/neurocluster/projects/biobank-cnv/BATCHES/CNVQC/referencefiles/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.bed");
		BufferedReader bedReader = new BufferedReader(new FileReader(CNVGoldStandard));
		
		dgvDels = new HashMap<String,IntervalTree<String>>();
		dgvDups = new HashMap<String,IntervalTree<String>>();

		String line;
		String data[];
		Pattern chrPattern = Pattern.compile("chr(\\S+)");
		int totalDels = 0;
		int totalDups = 0;
		
		while ((line = bedReader.readLine()) != null) {
			
			data = line.split("\t");
			
			String chr = data[0];
			Matcher chrMatcher = chrPattern.matcher(chr);
			if (chrMatcher.matches()) {
				chr = chrMatcher.group(1);
			}
			int start = Integer.valueOf(data[1]);
			int end = Integer.valueOf(data[2]);
			String id = data[3];
			String subType = data[6];
						
			if (subType.equals("Gain")) {
				
				totalDups++;
				if (dgvDups.containsKey(chr)) {
					dgvDups.get(chr).put(start, end, id);
				} else {
					IntervalTree<String> chrTree = new IntervalTree<String>();
					chrTree.put(start, end, id);
					dgvDups.put(chr, chrTree);
				}
				
			} else if (subType.equals("Loss")) {
				
				totalDels++;
				if (dgvDels.containsKey(chr)) {
					dgvDels.get(chr).put(start, end, id);
				} else {
					IntervalTree<String> chrTree = new IntervalTree<String>();
					chrTree.put(start, end, id);
					dgvDels.put(chr, chrTree);
				}
				
			}

		}
			
		bedReader.close();
		
		System.err.println("Total Gold Standard DELs loaded: " + totalDels);
		System.err.println("Total Gold Standard DUPs loaded: " + totalDups);
		
	}
	
	public List<CNV> getAllCNVs(int lineStart, int lineEnd) throws NumberFormatException, IOException {
		
		List<CNV> cnvs = new ArrayList<CNV>();
		String line;
		
		int totalCNVs = 0;
		int lineNum = 1;
		
		while ((line = cnvReader.readLine()) != null) {
			if (lineNum > lineEnd) {
				break;
			} else if (lineNum >= lineStart) {
				CNV currentCNV = parseCNVLine(line);
				if (currentCNV != null) {
					totalCNVs++;
					cnvs.add(currentCNV);
				}
			}
			lineNum++;
			
		}
		
		System.err.println("Total CNVs loaded: " + totalCNVs);
		return cnvs;
		
	}
	private CNV parseCNVLine(String line) throws NumberFormatException, IOException {
		
		//PennCNV is white-space delim... whereas my LRR BAF info isn't
		String lrrbaf[] = line.split("\t");
		String data[] = lrrbaf[0].split("\\s+");
		
		//New Code to parse CNV:
		//Only want individuals that passed AFFY QC
		String splitFile = data[4];

		if (sampleInformation.containsKey(splitFile)) {				
		
			SampleInformation sampInfo = sampleInformation.get(splitFile);
			
			String chr = null;
			int start = -1;
			int end = -1;
			Matcher locationParser = Pattern.compile("chr([\\dXY]{1,2}):(\\d+)\\-(\\d+)").matcher(data[0]);
			if (locationParser.matches()) {
				chr = locationParser.group(1);
				start = Integer.parseInt(locationParser.group(2));
				end = Integer.parseInt(locationParser.group(3));
			}
	
			int probeCount = (int)EqualSpliter(data[1]);
			int copyNumber = (int)EqualSpliter(data[3]);
			double conf = EqualSpliter(data[7]);
			//New code to parse CNV
						
			if (!samples.contains(splitFile)) {
				samples.add(splitFile);
			}
			if (!chrs.contains(chr)) {
				chrs.add(chr);
			}

			CNV currentCNV = new CNV(
					chr,
					start,
					end,
					copyNumber,
					probeCount,
					conf,
					sampInfo,
					getLRRBAF(sampInfo.getSplitFile(), chr, start, end)
					);
			
			currentCNV.setGoldStandardSV(GoldStandardParser(chr, start, end, currentCNV.getCopyType()));
			currentCNV.setDistTel(getDistTelo(chr, start, end));
			currentCNV.setDistCen(getDistCen(chr, start, end));
			
			// Checks for intersection to WES baits -- even if the CNV doesn't have WES data (for annotation purposes)
			IntervalTree<Boolean> wesBaits = getIntersectingBaits(chr, start, end);	
			currentCNV.setTotalIntersectingBaits(wesBaits.size());
			
			return currentCNV;
			
		} else {
			
			return null;
			
		}

	}
	private IntervalTree<Boolean> getIntersectingBaits(String chr, int start, int end) throws IOException {
		
		IntervalTree<Boolean> baitTree = new IntervalTree<Boolean>();
		
		Iterator itr = baitTabix.query(chr, start, end);

		String line;
		
		while ((line = itr.next()) != null) {

			String data[] = line.split("\t");
			
			int baitStart = Integer.parseInt(data[1]);
			int baitEnd = Integer.parseInt(data[2]);
			
			baitTree.put(baitStart, baitEnd, false);
			
		}
				
		return baitTree;
		
	}
	
	private LRRandBAFInformation getLRRBAF(File splitFile, String chr, int start, int end) throws IOException {
		
		TabixReader lrrbafTabixReader = new TabixReader(splitFile.getAbsolutePath(), splitFile.getAbsolutePath() + ".tbi");
		
		int len = end - start;
		int qStart = (start - len) < 0 ? 0 : (start - len);
		int qEnd = end + len;
	
		Iterator itr = lrrbafTabixReader.query(chr, qStart, qEnd);
		
		String line;
		String data[];
		int nLeft = 0;
		int nRight = 0;
				
		DescriptiveStatistics lrrStat = new DescriptiveStatistics();
		DescriptiveStatistics bafStat = new DescriptiveStatistics();
		
		while ((line = itr.next()) != null) {
			
			data = line.split("\t");
			int currPos = Integer.parseInt(data[1]);
			if (currPos < start) {
				nLeft++;
			} else if (currPos > (end + 1)) {
				nRight++;
			} else {
				double lrr = Double.parseDouble(data[4]);
				double baf = Double.parseDouble(data[5]);
		
				lrrStat.addValue(lrr);
				bafStat.addValue(Math.abs(0.5 - baf));
			}
			
		}
				
		LRRandBAFInformation lrrbaf = new LRRandBAFInformation(lrrStat, bafStat, nLeft, nRight);
		
		return lrrbaf;
		
	}
	public class LRRandBAFInformation {
		
		private double lrrMean;
		private double lrrSD;
		private double lrrMedian;
		private double bafMean;
		private double bafSD;
		private double bafMedian;
		private int nLeft;
		private int nRight;
		private DecimalFormat df;
		
		private LRRandBAFInformation(DescriptiveStatistics lrrStat, DescriptiveStatistics bafStat, int nLeft, int nRight) {
			lrrMean = lrrStat.getMean();
			lrrSD = lrrStat.getStandardDeviation();
			lrrMedian = lrrStat.getPercentile(0.50);
			bafMean = bafStat.getMean();
			bafSD = bafStat.getStandardDeviation();
			bafMedian = bafStat.getPercentile(0.50);
			this.nLeft = nLeft;
			this.nRight = nRight;
			df = new DecimalFormat("##.#######");
		}
		
		private String format(double val) {
			if (Double.isNaN(val)) {
				return "NaN";
			} else {
				return df.format(val);
			}
		}
		public String returnPrintable() {
			return (format(bafMean) + "\t" + format(bafSD) + "\t" + format(bafMedian) + "\t" + format(lrrMean) + "\t" + format(lrrSD) + "\t" + format(lrrMedian) + "\t" + nLeft + "\t" + nRight);
		}
		
	}
	
	public String GoldStandardParser (String chr, int start, int end, CopyType ct) {
		
		IntervalTree<String> toParse;
		if (ct == CopyType.DEL) {
			toParse = dgvDels.get(chr);
		} else {
			toParse = dgvDups.get(chr);
		}

		String maxHit = null;
		double maxVal = 0;
		
		Interval cnvInt = new Interval(chr, start, end);
		java.util.Iterator<Node<String>> hitCNVs = toParse.overlappers(start, end);
		
		while (hitCNVs.hasNext()) {
			Node<String> hitCNV = hitCNVs.next();
			Interval gsCnvInt = new Interval(chr, hitCNV.getStart(), hitCNV.getEnd());
			
			double compOne = gsCnvInt.getIntersectionLength(cnvInt);
			double overlapOne = compOne / (double) gsCnvInt.length();
			
			double compTwo = cnvInt.getIntersectionLength(gsCnvInt);
			double overlapTwo = compTwo / (double) cnvInt.length();
			
			double score = overlapOne * overlapTwo;
			
			if (overlapOne > 0.8 && overlapTwo > 0.8) {
			
				if (score > maxVal) {
					maxVal = score;
					maxHit = hitCNV.getValue();
				}
			}
		}
		
		return maxHit;
		
	}

	private Map<String, CytobandInfo> buildCytoband(File cytoband) throws IOException {
		
		Map<String, CytobandInfo> cytoMap = new HashMap<String, CytobandInfo>();
		
		BufferedReader cytoReader = new BufferedReader(new FileReader(cytoband));
		
		String line;
		String data[];
		
		while ((line = cytoReader.readLine()) != null) {
			
			data = line.split("\t");
			cytoMap.put(data[0], new CytobandInfo(Integer.parseInt(data[1]),Integer.parseInt(data[4]), Integer.parseInt(data[2]), Integer.parseInt(data[3])));
			
		}
		
		cytoReader.close();
		
		return cytoMap;
		
	}
	
	private double getDistTelo(String chr, int start, int end) {
		
		CytobandInfo cyto = cytoMap.get(chr);
		
		//Start by getting distance to the telomere from either end of the CNV
		double startDist = start - cyto.getChrStart();
		double endDist = cyto.getChrStop() - end;
		
		// Then calculate smallest amount
		double minDist = Math.min(startDist, endDist);
			
		// Then return the %tage position on the chromosome
		return minDist / (double) end;
				
	}
	private double getDistCen(String chr, int start, int end) {
		
		CytobandInfo cyto = cytoMap.get(chr);
		
		int leftCentDist = Math.min(Math.abs(start - cyto.getCentLeft()),Math.abs(end - cyto.getCentLeft()));
		int rightCentDist = Math.min(Math.abs(start - cyto.getCentRight()),Math.abs(end - cyto.getCentRight()));
		
		if (leftCentDist < rightCentDist) {
			return leftCentDist / (double) cyto.getCentLeft();
		} else {
			return rightCentDist / ( (double) cyto.getChrStop() - (double) cyto.getCentRight());
		}
				
	}
	
	private class CytobandInfo {
		
		private int chrStart;
		private int chrStop;
		private int centLeft;
		private int centRight;
		
		private CytobandInfo(int chrStart, int chrStop, int centLeft, int centRight) {
			this.chrStart = chrStart;
			this.chrStop = chrStop;
			this.centLeft = centLeft;
			this.centRight = centRight;
		}

		public int getChrStart() {
			return chrStart;
		}
		public int getChrStop() {
			return chrStop;
		}
		public int getCentLeft() {
			return centLeft;
		}
		public int getCentRight() {
			return centRight;
		}
		
	}
	
	//Utilities to parse raw CNVs
	private double EqualSpliter (String toParse) {
		String parsed[] = toParse.split("\\=");
		String replaced = parsed[1].replaceAll(",","");
		return Double.parseDouble(replaced);
	}
	
}
