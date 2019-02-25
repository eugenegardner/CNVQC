package merger;

import java.io.File;
import java.io.IOException;
import java.util.List;

import merger.CNVConverter.CopyType;
import merger.cluster.CNVClusterer;
import merger.cluster.ReciprocalClustering;
import sampleannotator.resources.SampleLoader;
import utilities.CNV;
import utilities.CNVAnnotator.OverlapError;
import utilities.CNVMergerOptions;

public class CNVMerger {
	
	private static File tmpDir;
	
	public CNVMerger(String args[]) throws IOException, OverlapError {
		
		CNVMergerOptions options = new CNVMergerOptions(args);
		File toMerge = options.getRawCNVs();
		tmpDir = options.getTmpDirectory();
		
		SampleLoader sampleLoader = new SampleLoader();
		
		//This will read raw CNVS and print all raw information necessary for filtering
		ProcessedCNVReader reader = new ProcessedCNVReader(toMerge, options.isFilter());
		List<CNV> CNVs = reader.getCNVs(sampleLoader.getSampleInformation()); //TRUE flag for only WES samples
		reader.close();
		
		//Build constructor for merging CNVs (sets up VCF writer and does initial CNV merge)
		
		CNVClusterer clusterer = new ReciprocalClustering(tmpDir);
		
		CNVMergerMethods mergerMethods = new CNVMergerMethods(options.getOutput(),
				options.getFastaRef(),
				options.getFastaIndex(),
				sampleLoader.getSamples(),
				reader.getChrs(), 
				tmpDir,
				CNVs,
				clusterer);
		
		mergerMethods.MergeCNVs(CopyType.DEL);
		mergerMethods.MergeCNVs(CopyType.DUP);
		
		mergerMethods.close();	
		
//		clusterer = new DBSCANClustering(options.getEpsilon(), tmpDir);
//		
//		mergerMethods = new CNVMergerMethods(options.getOutput(),
//				options.getFastaRef(),
//				options.getFastaIndex(), 
//				reader.getSamples(), 
//				reader.getChrs(), 
//				tmpDir,
//				CNVs,
//				clusterer);
//		
//		mergerMethods.MergeCNVs(CopyType.DEL);
//		mergerMethods.MergeCNVs(CopyType.DUP);
//		
//		mergerMethods.close();	
		
	}
	
}
