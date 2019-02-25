package loaders;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.IntervalTree;

public class PathogenicLoader {

	private Map<String, IntervalTree<PathogenicLocus>> pathogenicDelIntervals;
	private Map<String, IntervalTree<PathogenicLocus>> pathogenicDupIntervals;

	
	public PathogenicLoader() throws IOException {
		
		BufferedReader pathogenicReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/loaders/pathogenic_cnvs.bed"), "UTF-8"));
		BuildPathogenicTrees(pathogenicReader);
		pathogenicReader.close();
	
	}
	
	public Map<String, IntervalTree<PathogenicLocus>> getPathogenicDupIntervals () {
		return pathogenicDupIntervals;
	}
	public Map<String, IntervalTree<PathogenicLocus>> getPathogenicDelIntervals () {
		return pathogenicDelIntervals;
	}
	
	private void BuildPathogenicTrees(BufferedReader pathogenicReader) throws IOException {
		
		Map<String, IntervalTree<PathogenicLocus>> dupLoci = new HashMap<String, IntervalTree<PathogenicLocus>>();
		Map<String, IntervalTree<PathogenicLocus>> delLoci = new HashMap<String, IntervalTree<PathogenicLocus>>();
		String line;
		String data[];
		
		for (int x = 1; x <= 22; x++) {
			dupLoci.put(String.valueOf(x), new IntervalTree<PathogenicLocus>());
			delLoci.put(String.valueOf(x), new IntervalTree<PathogenicLocus>());
		}
		dupLoci.put("X", new IntervalTree<PathogenicLocus>());
		dupLoci.put("Y", new IntervalTree<PathogenicLocus>());
		delLoci.put("X", new IntervalTree<PathogenicLocus>());
		delLoci.put("Y", new IntervalTree<PathogenicLocus>());
		
		while ((line = pathogenicReader.readLine()) != null) {
			
			data = line.split("\t");
			String chr = data[0];
			int start = Integer.parseInt(data[1]);
			int stop = Integer.parseInt(data[2]);
			String name = data[3];
			int copyNumber = Integer.parseInt(data[4]);
			Set<String> genes = getGenes(data[6]);
			IntervalTree<Integer> exons = getExons(data[7], data[8], data[9]);
			PathogenicType pType = PathogenicType.valueOf(data[5]);
						
			int length = stop - start;
			PathogenicLocus locus = new PathogenicLocus(name, copyNumber, genes, pType, length, exons);
			if (copyNumber == 3) {
				dupLoci.get(chr).put(start, stop, locus);
			} else {
				delLoci.get(chr).put(start, stop, locus);
			}
						
		}

		pathogenicDelIntervals = new HashMap<String, IntervalTree<PathogenicLocus>>(delLoci);
		pathogenicDupIntervals = new HashMap<String, IntervalTree<PathogenicLocus>>(dupLoci);
		delLoci.clear();
		dupLoci.clear();
		
	}
		
	public class PathogenicLocus {
		
		private String name;
		private int copyNumber;
		private Set<String> requiredGenes;
		private PathogenicType pathogenicType;
		private int length;
		private IntervalTree<Integer> exons;
		
		public PathogenicLocus(String name, int copyNumber, Set<String> requiredGenes, PathogenicType pathogenicType, int length, IntervalTree<Integer> exons) {
			
			this.name = name;
			this.copyNumber = copyNumber;
			this.requiredGenes = requiredGenes;
			this.pathogenicType = pathogenicType;
			this.length = length;
			this.exons = exons;
						
		}

		public String getName() {
			return name;
		}
		public int getCopyNumber() {
			return copyNumber;
		}
		public Set<String> getRequiredGenes() {
			return requiredGenes;
		}
		public PathogenicType getPathogenicType() {
			return pathogenicType;
		}
		public int getLength () {
			return length;
		}
		public IntervalTree<Integer> getExons() {
			return exons;
		}
		
	}
	
	public enum PathogenicType {
		FIFTY,GENEREQUIRED,ONEHUNDRED,EXONS,DUP,MB,MBGENE,SEGDUP;
	}

	private Set<String> getGenes (String geneList) {
		Set<String> genes = new HashSet<String>();
		if (!geneList.equals("null")) {
			if (geneList.length() > 0) {
				String splitGenes[] = geneList.split(",");
				for (String g : splitGenes) {
					genes.add(g);
				}
			}
		}
		return genes;
	}
	private IntervalTree<Integer> getExons (String startString, String exonLengths, String exonStarts) {
		
		IntervalTree<Integer> exons = new IntervalTree<Integer>();

		if (!startString.equals("null")) {
		
			int geneStart = Integer.parseInt(startString);
			String[] lengthsArray = exonLengths.split(",");
			String[] startsArray = exonStarts.split(",");
			
			for (int x = 0; x <= (startsArray.length - 1); x++) {
					
				int currentStart = geneStart + Integer.parseInt(startsArray[x]);
				int currentStop = currentStart + Integer.parseInt(lengthsArray[x]);
					
				exons.put(currentStart, currentStop, x + 1);
					
			}
		}
		
		return exons;
	}
	
}
