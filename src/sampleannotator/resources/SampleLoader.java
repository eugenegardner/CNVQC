package sampleannotator.resources;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sampleannotator.resources.SampleInformation.Gender;

public class SampleLoader {

	private Map<String, SampleInformation> sampleInformation;
	private Set<String> samples;
	
	public SampleLoader() throws IOException {
		
		samples = new HashSet<String>();
		
		BufferedReader affyReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/sampleannotator/resources/affy_sample_table.csv"), "UTF-8"));
		//Do initially as affy2CEL to ensure we only get one CEL file per individual
		sampleInformation = BuildAffyCelRelationship(affyReader);
		affyReader.close();

		BufferedReader keyReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/sampleannotator/resources/sample_name_key.csv"), "UTF-8"));
		addWESInfo(keyReader);
		keyReader.close();
		
		//Switch to CEL2affy so we can access info from splitReader (and make initial SampleInfo hash).
		sampleInformation = InvertSampleID();
		
		BufferedReader splitReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream("/sampleannotator/resources/CELtoSplit.txt"), "UTF-8"));

		addSplitFileAndQCData(splitReader);
		
		splitReader.close();
		System.err.println("Total samples parsed for sample information: " + sampleInformation.size());
		
	}
	
	public Map<String, SampleInformation> getSampleInformation() {
		return sampleInformation;
	}
	public Set<String> getSamples() {
		return samples;
	}
	
 	private Map<String, SampleInformation> BuildAffyCelRelationship(BufferedReader affyReader) throws IOException {
		
		Map<String, SampleInformation> affy2CEL = new HashMap<String, SampleInformation>();
		Map<String, Double> qualityScores = new HashMap<String, Double>();
		
		String line;
		String data[];
		int dropped = 0;
		
		while ((line = affyReader.readLine()) != null) {
		
			data = line.split(",");
			
			Double quality = determineQuality(data[17], data[18], data[19]);
			String affyID = data[5];
			String CELID = data[9] + ".CEL";
			boolean pass = data[21].equals("Pass");
			boolean ctrl = affyID.contains("HG");
			
			String reportedGender = data[23].toUpperCase();
			String actualGender = data[24].toUpperCase();
			Gender gender = getGender(reportedGender, actualGender);
			
			if (pass && !ctrl) {
				if (qualityScores.containsKey(affyID)) {
					if (qualityScores.get(affyID) < quality) {
						affy2CEL.put(affyID, new SampleInformation(CELID, affyID, gender, quality));
						qualityScores.put(affyID, quality);
					} else {
						dropped++;
					}
				} else {
					affy2CEL.put(affyID, new SampleInformation(CELID, affyID, gender, quality));
					qualityScores.put(affyID, quality);
				}
			} else {
				dropped++;
			}
		}
		
		System.err.println("Samples dropped due to duplication or quality fail: " + dropped);

		return affy2CEL;
		
	}
 	
 	private double determineQuality(String dQCString, String QCCRString, String clusterCRString) {
 		
 		double dQC;
		try {
			dQC = Double.parseDouble(dQCString);
		} catch (NullPointerException | NumberFormatException e) {
			dQC = Double.NaN;
		}
 		
		double QCCR;
		try {
			QCCR = Double.parseDouble(QCCRString);
		} catch (NullPointerException | NumberFormatException e) {
			QCCR = Double.NaN;
		}
		
		double clusterCR;
		try {
			clusterCR = Double.parseDouble(clusterCRString);
		} catch (NullPointerException | NumberFormatException e) {
			clusterCR = Double.NaN;
		}
		
 		if (!Double.isNaN(clusterCR)) {
 			return clusterCR / 100; 
 		} else if (!Double.isNaN(QCCR)) {
 			return QCCR / 100;
 		} else if (!Double.isNaN(dQC)) {
 			return clusterCR;
 		} else {
 			return 0.0;
 		}
 		
 	}
 	
 	private Gender getGender (String reportedGender, String actualGender) {
 		Gender gender;
		if (reportedGender.equals(actualGender)) {
			gender = Gender.valueOf(actualGender);
		} else if (reportedGender.equals("Unknown")) {
			if (actualGender.equals("")) {
				gender = Gender.UNKNOWN;
			} else {
				gender = Gender.valueOf(actualGender);
			}
		} else if (actualGender.equals("")) {
			gender = Gender.valueOf(reportedGender);
		} else {
			gender = Gender.valueOf(actualGender);
		}
		return gender;
 	}

	private void addWESInfo(BufferedReader keyReader) throws IOException {
		
		String line;
		String data[];
		int WESsamples = 0;
				
		while ((line = keyReader.readLine()) != null) {
		
			data = line.split(",");
			String EGAN = data[2];
			String affyID = data[1];
			String sangerID = data[0];
			SampleInformation sampleInfo = sampleInformation.get(affyID);
			if (sampleInfo == null) {
				continue;
			}
			sampleInfo.setEGAN(EGAN);
			sampleInfo.setSangerID(sangerID);
			sampleInformation.put(affyID, sampleInfo);
			WESsamples++;
			
		}
		
		System.err.println("Total samples with WES information attached: " + WESsamples);
		
	}
	private Map<String, SampleInformation> InvertSampleID() {
		
		Map<String, SampleInformation> invertedSampleInformation = new HashMap<String, SampleInformation>();
		for (Map.Entry<String, SampleInformation> sampleEntry : sampleInformation.entrySet()) {
			invertedSampleInformation.put(sampleEntry.getValue().getCELName(), sampleEntry.getValue());
		}
		return invertedSampleInformation;
		
	}
	
	//This is currently hard coded. Might be bad!
	private Map<File, String[]> buildQCInformation() throws IOException {
		
		Map<File, String[]> qcData = new HashMap<File, String[]>();
		
		for (int x = 1; x <= 10; x++) {
			
			BufferedReader qcReader = new BufferedReader(new FileReader(new File("/lustre/scratch115/projects/interval_cnv/calling/batch" + x + "/batch" + x + ".qcsum")));
			
			String line;
			String data[];
			
			while ((line = qcReader.readLine()) != null) {
				data = line.split("\t");
				if (data[0].equals("File")) {
					continue;
				} else {
					File splitFile = new File(data[0]);
					qcData.put(splitFile, Arrays.copyOfRange(data, 1, data.length));
				}
			}
			
			qcReader.close();
						
		}
		
		return qcData;
		
	}
	//Call rate is current FUCKED -- Unsure if this is taken care of by QC file provided by Tao
	private void addSplitFileAndQCData(BufferedReader splitReader) throws IOException {
		
		String line;
		String data[];
		
		Map<String, SampleInformation> finalSampleInformation = new HashMap<String, SampleInformation>();
		
		Pattern filePattern = Pattern.compile("(split\\d+\\.a\\d{6}\\S*)");
		
		Map<File, String[]> qcInfo = buildQCInformation();
		
		while ((line = splitReader.readLine()) != null) {
			
			data = line.split("\t");
			String CELName = new File(data[0]).getName();
			File splitFile = new File(data[1]);
			if (sampleInformation.get(CELName) != null) {
 				
				SampleInformation sampInfo = sampleInformation.get(CELName);
 				sampInfo.setSplitFile(splitFile);
				String rawQCStats[] = qcInfo.get(splitFile);
				
				double wf = Double.parseDouble(rawQCStats[7]);
				int numCNV = Integer.parseInt(rawQCStats[8]);
				
				sampInfo.addFilterInformation(wf, numCNV);
				
				sampInfo.setLrr_mean(Double.parseDouble(rawQCStats[0]));
				sampInfo.setLrr_median(Double.parseDouble(rawQCStats[1]));
				sampInfo.setLrr_sd(Double.parseDouble(rawQCStats[2]));
				
				sampInfo.setBaf_mean(Double.parseDouble(rawQCStats[3]));
				sampInfo.setBaf_median(Double.parseDouble(rawQCStats[4]));
				sampInfo.setBaf_sd(Double.parseDouble(rawQCStats[5]));
				
				sampInfo.setBaf_drift(Double.parseDouble(rawQCStats[6]));
								
				finalSampleInformation.put(splitFile.getName(), sampInfo);
				
				// Add all samples that pass QC into sample list
				Matcher fileMatcher = filePattern.matcher(splitFile.getName());
				
				if (fileMatcher.matches()) {
					String fileName = fileMatcher.group(1);
					if (!samples.contains(fileName)) {
						samples.add(fileName);
					}
				}
				
			}
			
			
			
		}
		
		sampleInformation = finalSampleInformation;
		
	}

}
