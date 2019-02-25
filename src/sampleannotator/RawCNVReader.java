package sampleannotator;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
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
		
		File axiomProbes = new File("/lustre/scratch115/projects/interval_cnv/calling/reference/axiom_probes.bed.gz");
		axiomTabix = new TabixReader(axiomProbes.getAbsolutePath(),axiomProbes.getAbsolutePath() + ".tbi");
		
		File baits = new File("/lustre/scratch115/projects/interval_cnv/calling/sanger_exome_cnvs/exome_baits/bait_regions_unpadded.merge500.bed.gz");
		baitTabix = new TabixReader(baits.getAbsolutePath(),baits.getAbsolutePath() + ".tbi");
		
		cytoMap = buildCytoband(new File("/lustre/scratch115/projects/interval_cnv/calling/reference/cytoBand.txt"));
		
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
		WESTabix.CONVEX.close();
		WESTabix.XHMM.close();
		WESTabix.CLAMMS.close();
		WESTabix.CANOES.close();

	}
	
	private void BuildGoldStandard() throws IOException {
		File CNVGoldStandard = new File("/lustre/scratch115/projects/interval_cnv/calling/reference/DGV.GS.March2016.50percent.GainLossSep.Final.hg19.bed");
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
	
	public CNV getRandomCNV() throws IOException {
		
		Random randomGen = new Random();
		int testCNV = randomGen.nextInt(32695);
		int lineNum = 0;
		String line;
		CNV cnv = null;
		while ((line = cnvReader.readLine()) != null) {
			if (lineNum == testCNV) {
				cnv = parseCNVLine(line);
				break;
			}
			lineNum++;
		}
		if (cnv == null) {
			return null;
		} else {
			if (cnv.getCopyNumber() < 2) {
				return cnv;
			} else {
				return null;
			}
		}
		
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
		Pattern filePattern = Pattern.compile("(split\\d+\\.a\\d{6}\\S*)");
		
		//New Code to parse CNV:
		//Only want individuals that passed AFFY QC
		File splitFile = new File(data[4]);

		if (sampleInformation.containsKey(splitFile.getName())) {				
		
			SampleInformation sampInfo = sampleInformation.get(splitFile.getName());
			
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
			
			//Captures sample names
			Matcher fileMatcher = filePattern.matcher(splitFile.getName());
			
			if (fileMatcher.matches()) {
				String fileName = fileMatcher.group(1);
				if (!samples.contains(fileName)) {
					samples.add(fileName);
				}
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
			
			//Check if there are any intersecting CNVs for this individual if individual has WES data
			if (sampInfo.hasWES()) {
					
//				System.out.println(currentCNV.getLocationCoordinates() + " -- Total Baits: " + probeCount);
				
//				System.out.println("CONVEX");
				currentCNV.setIntersectingWESConvexCNVs(getWESCNVs(sampInfo.getEGAN(), currentCNV, WESTabix.CONVEX, wesBaits)); //This checks CONVEX CNVs.
//				System.out.println("XHMM");
				currentCNV.setIntersectingWESXHMMCNVs(getWESCNVs(sampInfo.getEGAN(), currentCNV, WESTabix.XHMM, wesBaits)); //This checks XHMM CNVs.
//				System.out.println("CLAMMS");
				currentCNV.setIntersectingWESCLAMMSCNVs(getWESCNVs(sampInfo.getSangerID(),currentCNV, WESTabix.CLAMMS, wesBaits)); //This checks CLAMMS CNVs.
//				System.out.println("CANOES");
				currentCNV.setIntersectingWESCANOESCNVs(getWESCNVs(sampInfo.getSangerID(),currentCNV, WESTabix.CANOES, wesBaits)); //This checks CLAMMS CNVs.
				
				if (currentCNV.getTotalIntersectingBaits() > 0) {
					currentCNV.setWESStats(getWESStats(chr, start, end, sampInfo.getEGAN()));
				}

			}

			return currentCNV;
			
		} else {
			
			return null;
			
		}

	}
	private IntervalTree<Boolean> getIntersectingProbes(String chr, int start, int end) throws IOException {
		
		IntervalTree<Boolean> baitTree = new IntervalTree<Boolean>();
		
		Iterator itr = axiomTabix.query(chr, start, end);

		String line;
		
		while ((line = itr.next()) != null) {

			String data[] = line.split("\t");
			
			int baitStart = Integer.parseInt(data[1]);
			int baitEnd = Integer.parseInt(data[2]);
			
			baitTree.put(baitStart, baitEnd, false);
			
		}
				
		return baitTree;
		
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
	
	private double getWESCNVs(String eganID, CNV cnv, WESTabix tabix, IntervalTree<Boolean> wesBaits) throws NumberFormatException, IOException {
		
		String line;
		String data[];
		
//		List<CNVInterval> WESCNVs = new ArrayList<CNVInterval>();
		
		Iterator itr = tabix.getReader().query(cnv.getChr(), cnv.getStart(), cnv.getEnd());			
		
		double totalBaits = wesBaits.size();
		double baitsIntersected = 0;
				
		double totalProbes = 0;
		double probesIntersected = 0;
								
		while ((line = itr.next()) != null) {
			
			data = line.split("\t");
				
			int startWES = Integer.parseInt(data[1]);
			int endWES = Integer.parseInt(data[2]);
			double numProbes = Double.parseDouble(data[3]);
			CopyType ctWES = CopyType.valueOf(data[5]);
			String eganIDWES = data[6];
			
			CNV wesCNV = new CNV(cnv.getChr(), startWES, endWES, 0, (int) numProbes, 0.0, null, null);
			ArrayList<CNV> toInt = new ArrayList<CNV>();
			toInt.add(cnv);
			toInt.add(wesCNV);
			Interval wesCnvInt = new Interval(cnv.getChr(), startWES, endWES);
			Interval cnvInt = new Interval(cnv.getChr(), cnv.getStart(), cnv.getEnd());
			
			// New trigger is if we have ±20% of probes overlapped by ArrayCNV
			// This means 80% of all probes that an array CNV overlaps should be within the WES interval AND the WES CNV contains no more than 20% extra probes

			if (eganIDWES.equals(eganID) && cnv.getCopyType().equals(ctWES) && wesCnvInt.intersects(cnvInt)) {
								
				// This checks for the WES CNV encompassing >80% of WES baits covered by the array CNV
				java.util.Iterator<Node<Boolean>> baitOverlap = wesBaits.overlappers(startWES, endWES);
								
				while (baitOverlap.hasNext()) {
					baitOverlap.next();
					baitsIntersected++;
				}
				
				// This checks for the Array CNV encompassing >80% of array probes covered by the WES CNV
				IntervalTree<Boolean> axiomProbes = getIntersectingProbes(cnv.getChr(), startWES, endWES);
				totalProbes += axiomProbes.size();		
				
				java.util.Iterator<Node<Boolean>> probeOverlap = axiomProbes.overlappers(cnv.getStart(), cnv.getEnd() + 1);
							
				while (probeOverlap.hasNext()) {
					probeOverlap.next();
					probesIntersected++;
				}
//								
//				double test1 = totalIntersectingBaits / totalBaits;
//				double test2 = totalIntersectingProbes / totalProbes;
//				System.out.println("\t -- Total WES Baits (" + totalBaits + "): " + totalIntersectingBaits + "\t(" + test1 + ")");
//				System.out.println("\t -- Total Axiom Pro (" + totalProbes + "): " + totalIntersectingProbes + "\t(" + test2 + ")");
//				
//				if (((test1) >= 0.80) && ((test2) >= 0.80)) {
//						
//						System.out.println("\tPASS");
//						WESCNVs.add(new CNVInterval(cnv.getChr(), startWES, endWES));
//				
//				}
			}
		}
		
//		System.out.println("\tWES: " + baitsIntersected + " -- " + wesBaits.size());
//		System.out.println("\tARR: " + probesIntersected + " -- " + totalProbes);
		
		if (totalProbes == 0) {
//			System.out.println(baitsIntersected / totalBaits);
			return(baitsIntersected / totalBaits);
		} else {
//			System.out.println((baitsIntersected / totalBaits) + (probesIntersected / totalProbes));
			return((baitsIntersected / totalBaits) * (probesIntersected / totalProbes));
		}
				
	}
	private DescriptiveStatistics getWESStats(String chr, int start, int end, String EGAN) throws IOException {
		
		File wesL2RFile = new File("/lustre/scratch115/projects/interval_cnv/calling/sanger_exome_cnvs/CONVEX/convex_out/ProbeRD_" + EGAN + "_LR2.coords.bed.gz");
		DescriptiveStatistics l2rStats = new DescriptiveStatistics();
		TabixReader l2rTabix = new TabixReader(wesL2RFile.getAbsolutePath(),wesL2RFile.getAbsolutePath() + ".tbi");
		
		Iterator itr = l2rTabix.query(chr, start, end);
		
		String line;
		String data[];
		
		while ((line = itr.next()) != null) {
			
			data = line.split("\t");
			double l2r = Double.parseDouble(data[3]);
			l2rStats.addValue(l2r);
			
		}
		return l2rStats;
		
	}
	
	private LRRandBAFInformation getLRRBAF(File splitFile, String chr, int start, int end) throws IOException {
		
		TabixReader lrrbafTabixReader = new TabixReader(splitFile.getAbsolutePath() + ".sorted.bed.gz", splitFile.getAbsolutePath() + ".sorted.bed.gz.tbi");
		
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
	
	private enum WESTabix {
		
		CONVEX(new File("/lustre/scratch115/projects/interval_cnv/calling/sanger_exome_cnvs/CONVEX/interval_cnv_calls.raw.ejg.sorted.bed.gz")),
		XHMM(new File("/lustre/scratch115/projects/interval_cnv/calling/sanger_exome_cnvs/XHMM/INTERVAL.xhmm.sorted.bed.gz")),
		CLAMMS(new File("/lustre/scratch115/projects/interval_cnv/calling/sanger_exome_cnvs/CLAMMS/CLAMMS.cnvs.bed.gz")),
		CANOES(new File("/lustre/scratch115/projects/interval_cnv/calling/sanger_exome_cnvs/CANOES/CANOES.cnvs.bed.gz"));
		
		private TabixReader reader;
		
		private WESTabix(File tabixFile) {
			try {
				reader = new TabixReader(tabixFile.getAbsolutePath(),tabixFile.getAbsolutePath() + ".tbi");
			} catch (IOException e) {
				System.err.println("Could not load WES tabix file: " + tabixFile.getAbsolutePath());
				e.printStackTrace();
				System.exit(1);
			}
		}

		public TabixReader getReader() {
			return reader;
		}
		public void close() {
			reader.close();
		}
				
	}
	
}
