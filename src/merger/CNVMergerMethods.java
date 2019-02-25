package merger;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import merger.CNVConverter.CopyType;
import merger.cluster.CNVClusterer;
import utilities.CNV;
import utilities.CNVAnnotator.OverlapError;
import utilities.CNVInterval;

public class CNVMergerMethods implements Closeable {

	private double epsilon;
	private VCFEngine vcfEngine;
	private File output;
	private File tmpDir;
	private List<CNV> cnvs;
	private CNVClusterer clusterer;
	private BufferedWriter perIndividualWriter;
	
	public CNVMergerMethods(File output, File fastaRef, File fastaIndex, Set<String> samples, Set<String> chromosomes, File tmpDir, List<CNV> cnvs, CNVClusterer clusterer) throws IOException {
		
		this.output = output;
		this.tmpDir = tmpDir;
		this.cnvs = cnvs;
		this.clusterer = clusterer;
		vcfEngine = new VCFEngine(new File(output.getAbsolutePath() + "." + clusterer.getClusterType()),
				fastaRef,
				fastaIndex,
				samples, 
				chromosomes);
		perIndividualWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + ".qcdMerged.txt")));
				
	}
	
	@Override
	public void close() throws IOException {
		vcfEngine.close();
		perIndividualWriter.close();
	}
	
	public List<CNV> MergeCNVs(CopyType ct) throws IOException, OverlapError {
		
		//Merge CNVs by simple overlap (just required > 0 overlap)
		BufferedWriter mergedWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + "." + ct + "." + clusterer.getClusterType() + ".merged.bed")));
		BufferedWriter rawWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + "." + ct + "." + clusterer.getClusterType() + ".raw_cnv.bed")));
		
		Map<String,IntervalTree<List<CNV>>> mergeIntervals = InitialMergerEngine.mergeCNVsbyOverlap(cnvs, ct, tmpDir);
		
		mergedWriter.write("track name=Merged" + ct + "-" + epsilon + " description=\"Merged CNVs INTERVAL\" itemRgb=\"On\"\n");
		rawWriter.write("track name=Original" + ct + "-" + epsilon + " description=\"Unmerged CNVs INTERVAL\" itemRgb=\"On\"\n");
		
		
		List<CNV> finalCNVs = new ArrayList<CNV>();
		
		int mergeNum = 1;
		
		for (Map.Entry<String, IntervalTree<List<CNV>>> currentChrEntry : mergeIntervals.entrySet()) {
			
			String chr = currentChrEntry.getKey();
			Iterator<Node<List<CNV>>> treeItr = currentChrEntry.getValue().iterator();
			
			int totalMerged = 0;
			while (treeItr.hasNext()) {

				Node<List<CNV>> currentNode = treeItr.next();
				
//				if (currentNode.getStart() == 34919625) {

				Map<Integer, CNVInterval> mergedCNVs = clusterer.getMergedCNVs(currentNode.getValue());
				
				for (Map.Entry<Integer, CNVInterval> cnvEntry : mergedCNVs.entrySet()) {
					
					CNVInterval currentInterval = cnvEntry.getValue();

					String color = checkColor(mergeNum);
										
					mergedWriter.write(chr + "\t" + currentInterval.getStart() + "\t" + currentInterval.getEnd() + "\t" + ct + "_" + mergeNum + "\t1000\t+\t" + currentInterval.getStart() + "\t" + currentInterval.getEnd() + "\t" + color + "\n");

					Map<String, IndividualRecord> individualCopyNumbers = new HashMap<String, IndividualRecord>();
					List<CNV> currentRawCNVs = currentInterval.getAttachedCNVs();
										
					for (CNV cnv : currentRawCNVs) {
					
						cnv.attachMergeGroup(mergeNum);
						finalCNVs.add(cnv);
						rawWriter.write(chr + "\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + ct + "_" + mergeNum + "-" + cnv.getSampleInformation().getSplitFile().getName() + "\t1000\t+\t" + cnv.getStart() + "\t" + cnv.getEnd() + "\t" + color + "\n");
						totalMerged++;
						individualCopyNumbers.put(cnv.getSampleInformation().getSplitFile().getName(), new IndividualRecord(cnv.getCopyNumber(),cnv.getConfidence()));
						perIndividualWriter.write(cnv.getPrintable() + "\t" + chr + ":" + currentInterval.getStart() + "-" + currentInterval.getEnd());
						perIndividualWriter.newLine();
						
					}
					
					perIndividualWriter.flush();
					vcfEngine.addRecord(chr, currentInterval.getStart(), currentInterval.getEnd(), individualCopyNumbers, ct);
					mergeNum++;
				}
				
				mergedWriter.flush();
				rawWriter.flush();
				
//			}
			}
			
			System.err.println("stats\t" + epsilon + "\t" + chr + "\t" + totalMerged);
			
		}
		
		mergedWriter.close();
		rawWriter.flush();
		rawWriter.close();
		
		return finalCNVs;
		
	}

	public class IndividualRecord {
		
		private int copyNumber;
		private double confidence;
				
		public IndividualRecord(int copyNumber, double confidence) {
			this.copyNumber = copyNumber;
			this.confidence = confidence;
		}

		public int getCopyNumber() {
			return copyNumber;
		}
		public double getConfidence() {
			return confidence;
		}
		
	}
	
	private String checkColor(int iteration) {
		
		Map<Integer, String> color = new HashMap<Integer, String>();
		color.put(0, "255,0,0");
		color.put(1, "255,127,0");
		color.put(2, "255,255,0");
		color.put(3, "0,255,0");
		color.put(4, "0,0,255");
		color.put(5, "75,0,130");
		color.put(6, "143,0,255");
		
		int mod = iteration % 7;
		return color.get(mod);
				
	}

	
}
