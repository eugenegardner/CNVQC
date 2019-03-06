package sampleannotator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;

import sampleannotator.RawCNVReader.LRRandBAFInformation;
import sampleannotator.resources.SampleInformation;
import sampleannotator.resources.SampleLoader;
import utilities.CNV;

public class CNVSampleAnnotator {

	private static BufferedWriter rawOutputWriter;
	
	public CNVSampleAnnotator(String args[]) throws IOException {
		
		CNVSampleAnnotatorOptions options = new CNVSampleAnnotatorOptions(args);
		File toAnnotate = options.getRawCNVs();
		
		SampleLoader sampleLoader = new SampleLoader(options.getQCSum(), options.getSplitDir());
		
		Map<String, SampleInformation> sampleInformation = sampleLoader.getSampleInformation();
		
		//This will read raw CNVS and print all raw information necessary for filtering
		RawCNVReader reader = new RawCNVReader(toAnnotate, sampleInformation);
		
		// Determine section of file to annotate;
		int fileStart;
		int fileEnd;
		
		int fileOption = options.getSectionOfFile();
		if (fileOption == -1) {
			//This means we are annotating the entire file
			fileStart = 1;
			fileEnd = Integer.MAX_VALUE;
		} else {
			int fileSection = fileOption - 1; // Have to do -1 to make the math work
			fileStart = 1 + (fileSection * 500);
			fileEnd = 500 + (fileSection * 500);
		}
		List<CNV> rawCNVs = reader.getAllCNVs(fileStart, fileEnd); //TRUE flag for only WES samples
		
		rawOutputWriter = new BufferedWriter(new FileWriter(new File(options.getOutput().getAbsolutePath() + ".txt")));
		PrintRawCNVs(rawCNVs);
		rawOutputWriter.close();
		reader.close();
		
	}
	
	private static void PrintRawCNVs(List<CNV> cnvs) throws IOException {
		
		DecimalFormat df = new DecimalFormat("##.##");
//		GenerateValidPlots validPlot = new GenerateValidPlots(new File("/lustre/scratch115/projects/interval_cnv/bed_files.lst"), tmpDir);
		
		printtab("#chr");
		printtab("start");
		printtab("end");
		printtab("SangerID");
		printtab("EGANID");
		printtab("split.file");
		printtab("has.wes");
		printtab("location");
		printtab("Copy_Number");
		printtab("Length_bp");
		printtab("Max_Log_BF");
		printtab("LRR_mean");
		printtab("LRR_median");
		printtab("LRR_SD");
		printtab("BAF_mean");
		printtab("BAF_median");
		printtab("BAF_SD");
		printtab("WF");
		printtab("BAF_drift");
		printtab("cel.file");
		printtab("NumCNV");
		printtab("Gender");
		printtab("density");
		printtab("callrate");
		printtab("wes.convex.int");
		printtab("wes.xhmm.int");
		printtab("wes.clamms.int");
		printtab("wes.canoes.int");
		printtab("wes.probe.count");
		printtab("DGVIntersect");
		printtab("WES.l2r.mean");
		printtab("WES.l2r.sd");
		printtab("num.l2r.probes");
		printtab("indiv.filtered");
		printtab("site.filtered");
		printtab("No_Probes");
		printtab("site.BAF_mean");
		printtab("site.BAF_SD");
		printtab("site.BAF_median");
		printtab("site.LRR_mean");
		printtab("site.LRR_SD");
		printtab("site.LRR_median");
		printtab("nLeft");
		printtab("nRight");
		printtab("abs.tel");
		printnewline("abs.cen");
				
		for (CNV cnv : cnvs) {

//			validPlot.GenerateCNVPlot(cnv);
			
			SampleInformation si = cnv.getSampleInformation();
			LRRandBAFInformation lrrbaf = cnv.getLRRBAF();
			
			printtab(cnv.getChr()); //0
			printtab(cnv.getStart()); //1
			printtab(cnv.getEnd()); //2
			printtab("NA"); //3
			printtab("NA"); //4
			printtab(si.getSplitFile().getAbsolutePath()); //5
			printtab("NA"); //6
			printtab(cnv.getLocationCoordinates()); //7
			printtab(cnv.getCopyNumber()); //8
			printtab(cnv.getLength()); //9
			printtab(df.format(cnv.getConfidence())); //10
			printtab(si.getLrr_mean()); //11
			printtab(si.getLrr_median()); //12
			printtab(si.getLrr_sd()); //13
			printtab(si.getBaf_mean()); //14
			printtab(si.getBaf_median()); //15
			printtab(si.getBaf_sd()); //16
			printtab(si.getWf()); //17
			printtab(si.getBaf_drift()); //18
			printtab("NA"); //19
			printtab(si.getNumCNV()); //20
			printtab(si.getGender()); //21
			printtab(df.format(cnv.getDensity())); //22
			printtab(df.format(si.getCallRate())); //23
			printtab("NA"); //24
			printtab("NA"); //25
			printtab("NA"); //26
			printtab("NA"); //27
			printtab(cnv.getTotalIntersectingBaits()); //28
			printtab(cnv.getGoldStandardSV()); //29
			printtab("NA"); //30
			printtab("NA"); //31
			printtab("NA"); //32
			printtab(si.isIndivFiltered()); //33
			printtab(cnv.isSiteFiltered()); //34
			printtab(cnv.getProbeCount()); //35
			printtab(lrrbaf.returnPrintable()); //36-43
			printtab(cnv.getDistTel());
			printnewline(cnv.getDistCen());
				
			rawOutputWriter.flush();

		}
		
	}
	public static void printtab(Object toPrint) throws IOException {
		rawOutputWriter.write(toPrint + "\t");
	}
	public static void printnewline(Object toPrint) throws IOException {
		rawOutputWriter.write(toPrint + "\n");
	}
	
}
