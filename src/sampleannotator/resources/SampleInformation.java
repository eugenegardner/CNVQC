package sampleannotator.resources;

import java.io.File;

public class SampleInformation {
	
	private File splitFile;
	private Gender gender;
	private Double callRate;
	private double lrr_sd;
	private double lrr_mean;
	private double lrr_median;
	private double baf_sd;
	private double baf_mean;
	private double baf_median;
	private double baf_drift;
	private double wf;
	private int numCNV;
	private boolean isIndivFiltered;
	
	public SampleInformation(Gender gender, 
			double callRate, 
			double lrr_sd, 
			double lrr_mean, 
			double lrr_median, 
			double baf_sd, 
			double baf_mean, 
			double baf_median, 
			double baf_drift, 
			double wf, 
			int numCNV,
			File splitFile) {
		
		this.gender = gender;
		this.callRate = callRate;
		this.lrr_sd = lrr_sd;
		this.lrr_mean = lrr_mean;
		this.lrr_median = lrr_median;
		this.baf_sd = baf_sd;
		this.baf_mean = baf_mean;
		this.baf_median = baf_median;
		this.baf_drift = baf_drift;
		this.wf = wf;
		this.numCNV = numCNV;
		this.splitFile = splitFile;
						
		isIndivFiltered = checkIndivFilter(numCNV, wf);

	}

	public File getSplitFile() {
		return splitFile;
	}
	public Gender getGender() {
		return gender;
	}
	public Double getCallRate() {
		return callRate;
	}
	public double getLrr_sd() {
		return lrr_sd;
	}
	public double getLrr_mean() {
		return lrr_mean;
	}
	public double getLrr_median() {
		return lrr_median;
	}
	public double getBaf_sd() {
		return baf_sd;
	}
	public double getBaf_mean() {
		return baf_mean;
	}
	public double getBaf_median() {
		return baf_median;
	}
	public double getBaf_drift() {
		return baf_drift;
	}
	public double getWf() {
		return wf;
	}
	public int getNumCNV() {
		return numCNV;
	}
	public boolean isIndivFiltered() {
		return isIndivFiltered;
	}
	
	private boolean checkIndivFilter (int numCNV, double wf) {
		
		if (numCNV <= 30 && wf > -0.03 && wf < 0.03) {
			return false;
		} else {
			return true;
		}
		
	}
	
	public enum Gender {
		MALE,
		FEMALE,
		UNKNOWN;		
	}


	
}


