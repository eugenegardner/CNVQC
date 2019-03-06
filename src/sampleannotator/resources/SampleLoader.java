package sampleannotator.resources;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import sampleannotator.resources.SampleInformation.Gender;

public class SampleLoader {

	private Map<String, SampleInformation> sampleInformation;
	private Set<String> samples;
	private File splitDir;
	
	public SampleLoader(File qcsum, File splitDir) throws IOException {
		
		samples = new HashSet<String>();
		
		this.splitDir = splitDir;
		
		BufferedReader splitReader = new BufferedReader(new FileReader(qcsum));
		sampleInformation = addSplitFileAndQCData(splitReader);
		
		splitReader.close();
		System.err.println("Total samples parsed for sample information: " + sampleInformation.size());
		
	}
	
	public Map<String, SampleInformation> getSampleInformation() {
		return sampleInformation;
	}
	public Set<String> getSamples() {
		return samples;
	}
	
	//Call rate is current FUCKED -- Unsure if this is taken care of by QC file provided by Tao
	private Map<String, SampleInformation> addSplitFileAndQCData(BufferedReader splitReader) throws IOException {
		
		String line;
		String data[];
		
		Map<String, SampleInformation> finalSampleInformation = new HashMap<String, SampleInformation>();
						
		while ((line = splitReader.readLine()) != null) {
			
			data = line.split("\t");
			
			String ID = data[0];
			File splitFile = new File(splitDir.getAbsolutePath() + "/" + ID + ".bed.gz");
						
			double lrr_mean = Double.parseDouble(data[1]);
			double lrr_median = Double.parseDouble(data[2]);
			double lrr_sd = Double.parseDouble(data[3]);
			
			double baf_mean = Double.parseDouble(data[4]);
			double baf_median = Double.parseDouble(data[5]);
			double baf_sd = Double.parseDouble(data[6]);
			
			double baf_drift = Double.parseDouble(data[7]);
							
			double wf = Double.parseDouble(data[8]);
			int numCNV = Integer.parseInt(data[9]);
			
			Gender gender = Gender.valueOf(data[10].toUpperCase());
			
			double callRate = Double.parseDouble(data[11]);

			if (callRate > 96.0) {
			
				finalSampleInformation.put(ID, new SampleInformation(gender, callRate, lrr_sd, lrr_mean, lrr_median, baf_sd, baf_mean, baf_median, baf_drift, wf, numCNV, splitFile));
			
				if (!samples.contains(ID)) {
					samples.add(ID);
				}
			}
			
		}
		
		return finalSampleInformation;
		
	}

}
