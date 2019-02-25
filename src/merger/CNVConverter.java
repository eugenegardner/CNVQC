package merger;

import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import utilities.CNV;
import utilities.CNVInterval;

public class CNVConverter {

	public CNVConverter() {
		super();
	}
	
	public static RealMatrix ListToMatrix(List<CNV> CNVs, boolean invert) {
		
		RealMatrix CNVMatrix = new Array2DRowRealMatrix(CNVs.size(),CNVs.size());
		for (int x = 0; x < CNVs.size(); x++) {
			
			CNV cnv1 = CNVs.get(x);
			int cnv1Start = cnv1.getStart();
			int cnv1End = cnv1.getEnd();
			CNVInterval cnv1Interval = new CNVInterval(cnv1.getChr(), cnv1Start, cnv1End);
			
			for (int y = 0; y < CNVs.size(); y++) {
			
				CNV cnv2 = CNVs.get(y);
				int cnv2Start = cnv2.getStart();
				int cnv2End = cnv2.getEnd();
				CNVInterval cnv2Interval = new CNVInterval(cnv2.getChr(), cnv2Start, cnv2End);
				
				double cnv1Overlap;
				double cnv2Overlap;
				
				if (cnv1Interval.intersects(cnv2Interval)) {
					double intersectLength = cnv1Interval.getIntersectionLength(cnv2Interval);
					cnv1Overlap = intersectLength / cnv1Interval.length();
					cnv2Overlap = intersectLength / cnv2Interval.length();
				} else {
					cnv1Overlap = 0.0;
					cnv2Overlap = 0.0;
				}
				
				if (invert) {
					cnv1Overlap = 1 - cnv1Overlap;
					cnv2Overlap = 1 - cnv2Overlap;
					CNVMatrix.setEntry(x, y, cnv1Overlap);
					CNVMatrix.setEntry(y, x, cnv2Overlap);
				} else {
					CNVMatrix.setEntry(x, y, cnv1Overlap);
					CNVMatrix.setEntry(y, x, cnv2Overlap);
				}
			}
		}
		
		return CNVMatrix;
		
	}
	
	public enum CopyType {
		DEL,
		DUP,
		OTHER;
	}
	
	public static CNVInterval GetCNVInterval(Set<CNV> currentMerge) {
		
		CNVInterval finalInterval = null;
		
		for (CNV cnv : currentMerge) {
			
			CNVInterval currentInterval = new CNVInterval(cnv.getChr(), cnv.getStart(), cnv.getEnd());
			
			if (finalInterval == null) {
				finalInterval = currentInterval;
			} else {
				finalInterval = finalInterval.union(currentInterval);
			}
			
		}
		
		return finalInterval;
	}
	
	//This is done as '<' because the matrix calculates similarity, NOT difference!
	public static boolean doubletonOverlap(List<CNV> CNVs) {
		
		RealMatrix matrix = CNVConverter.ListToMatrix(CNVs,true);
		double overlap1 = matrix.getEntry(0, 1);
		double overlap2 = matrix.getEntry(1, 0);
		return (overlap1 < 0.5 && overlap2 < 0.5);
		
	}
	
}
