package merger.cluster;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.ExecuteException;
import org.apache.commons.math3.linear.RealMatrix;

import merger.CNVConverter;
import utilities.CNV;
import utilities.CNVInterval;
import utilities.Combine;

public class DBSCANClustering implements CNVClusterer {
	
	private String chr; // Set during runtime for error reporting, nothing else.
	private int pos; // Set during runtime for error reporting, nothing else.
	private double eps;
	private File tmpDir;
	
	private final ClusterType clusterType = ClusterType.DBSCAN;
	
	public DBSCANClustering (double eps, File tmpDir) throws IOException {
		
		this.eps = eps;
		this.tmpDir = tmpDir;
		
	}
	
	public Map<Integer, CNVInterval> getMergedCNVs(List<CNV> currentCNVList) throws IOException {
		
		Map<Integer, CNVInterval> finalCNVs;
		
		if (currentCNVList.size() == 1) {
			finalCNVs = new HashMap<Integer, CNVInterval>();
			CNV currentCNV = currentCNVList.get(0);
			finalCNVs.put(1, new CNVInterval(currentCNV.getChr(), currentCNV.getStart(), currentCNV.getEnd(), currentCNV));
		} else if (currentCNVList.size() == 2) {
			
			finalCNVs = new HashMap<Integer, CNVInterval>();
			
			CNVInterval intOne = new CNVInterval(currentCNVList.get(0).getChr(), currentCNVList.get(0).getStart(), currentCNVList.get(0).getEnd(), currentCNVList.get(0));
			CNVInterval intTwo = new CNVInterval(currentCNVList.get(1).getChr(), currentCNVList.get(1).getStart(), currentCNVList.get(1).getEnd(), currentCNVList.get(1));
						
			if (CNVConverter.doubletonOverlap(currentCNVList)) {
				CNVInterval finalInterval = intOne.union(intTwo);
				finalCNVs.put(1, finalInterval);
			} else {
				finalCNVs.put(1, intOne);
				finalCNVs.put(2, intTwo);
				
			}
			
		} else {
			Map<String, CNV> CNVsToMerge = PrintMatrix(currentCNVList, CNVConverter.ListToMatrix(currentCNVList, true));
			runR();
			finalCNVs = parseResults(CNVsToMerge);
		}
		
		return finalCNVs;
	}
	public ClusterType getClusterType() {
		return clusterType;
	}
	
 	private Map<Integer, CNVInterval> parseResults(Map<String, CNV> CNVMap) throws IOException {
		BufferedReader resultsReader = new BufferedReader(new FileReader(new File(tmpDir.getAbsolutePath() + "/DBSCANClust.results.txt")));
		String line;
		String data[];
		
		Map<Integer, CNVInterval> mergedCNVs = new HashMap<Integer, CNVInterval>();
		List<CNV> soloMergedCNVs = new ArrayList<CNV>();
		
		while ((line = resultsReader.readLine()) != null) {
			data = line.split("\t");
			CNV currCNV = CNVMap.get(data[0]);

			CNVInterval currentInterval = new CNVInterval(data[0], currCNV.getStart(), currCNV.getEnd(), currCNV);
			int cluster = Integer.parseInt(data[1]);

			if (cluster == 0) {
				soloMergedCNVs.add(currCNV);
			} else {
				if (mergedCNVs.containsKey(cluster)) {
					CNVInterval adjustedInterval = mergedCNVs.get(cluster);
					CNVInterval unionInterval = adjustedInterval.union(currentInterval);
					mergedCNVs.put(cluster, unionInterval);
				} else {
					mergedCNVs.put(cluster, currentInterval);
				}
			}
		}
		resultsReader.close();
		CleanRData();
		return mergeLists(mergedCNVs, soloMergedCNVs);
				
	}
	private void runR() throws IOException {
		
		String osName = System.getProperty("os.name");
		//First cmd-line is OSX, second is *nix
		String command;
		if (osName.contains("OS X")) {
			command = "/Library/Frameworks/R.framework/Resources/bin/Rscript --vanilla /Users/eg15/Documents/Current Projects/INTERVAL/RawCNVCalls/merge_Java.R " + tmpDir.getAbsolutePath() + "/DBSCANClust.mat " + tmpDir.getAbsolutePath() + "/DBSCANClust.results.txt " +  eps;
		} else { //anything else?
			command = "/nfs/users/nfs_e/eg15/EugeneTools/R-3.5.1/bin/Rscript --vanilla /nfs/ddd0/eg15/merge_Java.R " + tmpDir.getAbsolutePath() + "/DBSCANClust.mat " + tmpDir.getAbsolutePath() + "/DBSCANClust.results.txt " +  eps;
		}
			
		CommandLine cmdLine = CommandLine.parse(command);
		DefaultExecutor executor = new DefaultExecutor();
		try {
			executor.execute(cmdLine);
		} catch (ExecuteException e) {
			e.printStackTrace();
			System.err.println("Offending merge: " + chr + "_" + pos);
			System.exit(1);
		}
		
	}
	private Map<String, CNV> PrintMatrix(List<CNV> CNVs, RealMatrix matrix) throws IOException {
		
		Map<String, CNV> mergedCNVs = new HashMap<String, CNV>();
		
		DecimalFormat df = new DecimalFormat("##.######");
		
		BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(tmpDir.getAbsolutePath() + "/DBSCANClust.mat")));
		
		//Header
		List<String> printList = new ArrayList<String>();
		printList.add("ID");
		for (int y = 0; y < matrix.getColumnDimension(); y++) {
			printList.add(CNVs.get(y).getSampleInformation().getSplitFile().getName() + "_" + CNVs.get(y).getChr() + "_" + CNVs.get(y).getStart() + "_" + CNVs.get(y).getEnd() + "_" + CNVs.get(y).getCopyNumber());
			mergedCNVs.put(CNVs.get(y).getSampleInformation().getSplitFile().getName() + "_" + CNVs.get(y).getChr() + ":" + CNVs.get(y).getStart() + "_" + CNVs.get(y).getEnd(), CNVs.get(y));
			chr = CNVs.get(y).getChr();
			pos = CNVs.get(y).getStart();
		}

		outWriter.write(Combine.combineList(printList, "\t") + "\n");
		//Actual Data
		for (int x = 0; x < matrix.getRowDimension(); x++) {
		
			printList = new ArrayList<String>();
			printList.add(CNVs.get(x).getSampleInformation().getSplitFile().getName() + "_" + CNVs.get(x).getChr() + "_" + CNVs.get(x).getStart() + "_" + CNVs.get(x).getEnd() + "_" + CNVs.get(x).getCopyNumber());
						
			for (int y = 0; y < matrix.getColumnDimension(); y++) {
				
				printList.add(df.format(matrix.getEntry(x, y)));
							
			}
			
			outWriter.write(Combine.combineList(printList, "\t") + "\n");
			
		}
		outWriter.close();
		return mergedCNVs;
		
	}
	private void CleanRData() {
		
		new File(tmpDir.getAbsolutePath() + "/DBSCANClust.mat").delete();
		new File(tmpDir.getAbsolutePath() + "/DBSCANClust.results.txt").delete();
		
	}
	
	private Map<Integer, CNVInterval> mergeLists(Map<Integer, CNVInterval> mergedCNVs, List<CNV> soloMergedCNVs) {
		
		Set<Integer> keys = mergedCNVs.keySet();
		int maxValue = keys.size() == 0 ? 0 : Collections.max(keys);
		
		for (CNV cnv : soloMergedCNVs) {
			maxValue++;
			mergedCNVs.put(maxValue, new CNVInterval(cnv.getChr(), cnv.getStart(), cnv.getEnd(), cnv));
		}
		
		return mergedCNVs;
	}

//	private List<CNVClusterable> PerformMDS(List<CompleteCNV> CNVs) throws IOException {
//		
//		RealMatrix matrix = CNVConverter.ListToMatrix(CNVs,true);
//		
//		double[][] input = matrix.getData();
//		double[][] output = MDSJ.stressMinimization(input);
//		int inputLength = input[0].length;
//		List<CNVClusterable> toCluster = new ArrayList<CNVClusterable>();
//		
//		BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + pos + "_" + ct + ".mds")));
//		outWriter.write("CNV\tX\tY\n");
//		
//		for (int x = 0; x < inputLength; x++) {
//			CompleteCNV currCNV = CNVs.get(x);
//			toCluster.add(new CNVClusterable(output[0][x], output[1][x], currCNV));
//			outWriter.write(currCNV.getSplitFile().getName() + "_" + currCNV.getChr() + "_" + currCNV.getStart() + "_" + currCNV.getEnd() + "\t" + output[0][x] + "\t" + output[1][x]);
//			outWriter.newLine();
//		}
//		
//		outWriter.close();
//		return toCluster;
//		
//	}
//	private void PerformDBSCAN(List<CNVClusterable> toCluster) throws IOException {
//		
//		DBSCANClusterer<CNVClusterable> dbscan = new DBSCANClusterer<CNVClusterable>(0.01, 2);
//		List<Cluster<CNVClusterable>> cluster = dbscan.cluster(toCluster);
//
//		BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(tmpDir.getAbsolutePath() + "/" + chr + "_" + pos + "_" + ct + ".clust")));
//		outWriter.write("CNV\tcluster\tX\tY\n");
//		
//		int clusterNum = 1;
//		int totalClustered = 0;
//		System.out.println("Num Clusters: " + cluster.size());
//				
//		for (Cluster<CNVClusterable> c : cluster) {
//			System.out.println("Cluster " + clusterNum + " size: " + c.getPoints().size());
//			totalClustered+=c.getPoints().size();
//			List<CNVClusterable> currCluster = c.getPoints();
//			for(CNVClusterable cnv : currCluster) {
//				
//				toCluster.remove(cnv);
//				CompleteCNV currCNV = cnv.getCNV();
//				CNVInterval currentInterval = new CNVInterval(currCNV.getChr(), currCNV.getStart(), currCNV.getEnd(), currCNV);
//				
//				if (mergedCNVs.containsKey(clusterNum)) {
//					CNVInterval adjustedInterval = mergedCNVs.get(clusterNum);
//					CNVInterval unionInterval = adjustedInterval.union(currentInterval);
//					mergedCNVs.put(clusterNum, unionInterval);
//					
//					List<CNVInterval> collection = mergedRawCNVs.get(clusterNum);
//					collection.add(currentInterval);
//					mergedRawCNVs.put(clusterNum, collection);
//				} else {
//					mergedCNVs.put(clusterNum, currentInterval);
//					List<CNVInterval> collection = new ArrayList<CNVInterval>();
//					collection.add(currentInterval);
//					mergedRawCNVs.put(clusterNum, collection);
//				}
//								
//				outWriter.write(currCNV.getSplitFile().getName() + "_" + currCNV.getChr() + "_" + currCNV.getStart() + "_" + currCNV.getEnd() + "\t" + clusterNum + "\t" + cnv.getPoint()[0] + "\t" + cnv.getPoint()[1]);
//				outWriter.newLine();
//			
//			}
//				
//			clusterNum++;
//
//		}
//				
//		System.out.println("Num Clustered: " + totalClustered);
//		
//				
//		for (CNVClusterable cnv : toCluster) {
//			CompleteCNV currCNV = cnv.getCNV();
//			
//			CNVInterval currentInterval = new CNVInterval(currCNV.getChr(), currCNV.getStart(), currCNV.getEnd(), currCNV);
//			soloMergedCNVs.add(currentInterval);
//						
//			outWriter.write(currCNV.getSplitFile().getName() + "_" + currCNV.getChr() + "_" + currCNV.getStart() + "_" + currCNV.getEnd() + "\t" + clusterNum + "\t" + cnv.getPoint()[0] + "\t" + cnv.getPoint()[1]);
//			outWriter.newLine();
//			clusterNum++;
//		}
//		System.out.println("Num Unclustered: " + toCluster.size());
//		outWriter.close();
//		
//	}
	
}
