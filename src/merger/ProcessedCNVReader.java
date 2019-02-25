package merger;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import sampleannotator.resources.SampleInformation;
import utilities.CNV;

public class ProcessedCNVReader implements Closeable {

	private BufferedReader processedReader;
	private Set<String> chrs;
	private boolean isFilter;
	
	public ProcessedCNVReader (File processedCNVs, boolean isFilter) throws IOException {
		
		processedReader = new BufferedReader(new FileReader(processedCNVs));
		chrs = new HashSet<String>();
		this.isFilter = isFilter;
		
	}
	
	public List<CNV> getCNVs(Map<String, SampleInformation> sampleInfo) throws IOException {
		
		String line;
		String data[];
		
		List<CNV> CNVs = new ArrayList<CNV>();
				
		while ((line = processedReader.readLine()) != null) {
			
			data = line.split("\t");
			
			String chr = data[0];
			int start = Integer.parseInt(data[1]);
			int end = Integer.parseInt(data[2]);
			File splitFile = new File(data[5]);
			int copyNumber = Integer.parseInt(data[8]);
			int probeCount = Integer.parseInt(data[35]);
			
			boolean passStatus = Boolean.parseBoolean(data[45]);
			double passValue = Double.parseDouble(data[60]);
						
			if (!chrs.contains(chr)) {
				chrs.add(chr);
			}

			SampleInformation sampleInformation = sampleInfo.get(splitFile.getName());
			if (isFilter) {
				if (passStatus) {
					CNVs.add(new CNV(chr, start, end, copyNumber, probeCount, passValue, sampleInformation, null, line));
				}
			} else {
				CNVs.add(new CNV(chr, start, end, copyNumber, probeCount, passValue, sampleInformation, null, line));
			}			
			
		}
		
		return CNVs;
		
	}
	public Set<String> getChrs() {
		return chrs;
	}
	
	@Override
	public void close() throws IOException {
		processedReader.close();
	}
	
}
