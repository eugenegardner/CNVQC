package merger;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import merger.CNVConverter.CopyType;
import merger.CNVMergerMethods.IndividualRecord;
import utilities.Combine;
import utilities.CNVAnnotator;
import utilities.CNVAnnotator.OverlapError;

public class VCFEngine implements Closeable {

	private BufferedWriter vcfWriter;
	private CNVAnnotator rawCNVAnnotator;
	private IndexedFastaSequenceFile fastaRef;
	private Set<String> samples;
	private Set<String> chrs;
	private DecimalFormat df;
	
	public VCFEngine(File output, File fastaRef, File fastaIndex, Set<String> samples, Set<String> chrs) throws IOException {
		
		vcfWriter = new BufferedWriter(new FileWriter(new File(output.getAbsolutePath() + ".vcf")));
		rawCNVAnnotator = new CNVAnnotator();
		this.fastaRef = new IndexedFastaSequenceFile(fastaRef, new FastaSequenceIndex(fastaIndex));
		this.samples = samples;
		this.chrs = chrs;
		df = new DecimalFormat("#.##");
		writeHeader();
		
	}
	
	@Override
	public void close() throws IOException {
		vcfWriter.close();
	}
	
	public void writeHeader() throws IOException {
		
		Date currentDate = GregorianCalendar.getInstance().getTime();
		DateFormat df = DateFormat.getDateTimeInstance();
		vcfWriter.write("##fileformat=VCFv4.2\n");
		vcfWriter.write("##fileDate=" + df.format(currentDate) + "\n");
		vcfWriter.write("##source=CNVPolisherv1.0\n");
		for (String chr : chrs) {
			int len = fastaRef.getSequence(chr).length();
			vcfWriter.write("##contig=<ID=" + chr + ",length=" + len + ">\n");
		}
		vcfWriter.write("##ALT=<ID=CNV,Description=\"Copy number variable region\">\n");
		vcfWriter.write("##ALT=<ID=DUP,Description=\"Duplication\">\n");
		vcfWriter.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles; If unknown, will be -1\">\n");
		vcfWriter.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
		vcfWriter.write("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count for this record\">\n");
		vcfWriter.write("##INFO=<ID=GENE,Number=.,Type=String,Description=\"Gene annotation for this variant\">\n");
		vcfWriter.write("##INFO=<ID=PATH,Number=.,Type=String,Description=\"Pathogenic annotation for this variant\">\n");
		vcfWriter.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		vcfWriter.write("##FORMAT=<ID=WC,Number=1,Type=Float,Description=\"WES Overlap Confidence by Random Forest\">\n");
		vcfWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + Combine.combineSet(samples, "\t") + "\n");
		vcfWriter.flush();
		
	}
	
	public void addRecord(String chr, int start, int end, Map<String, IndividualRecord> cnvs, CopyType ct) throws IOException, OverlapError {
		
		List<String> genotypeList = new ArrayList<String>();

		int ac = 0;
		
		for (String sample : samples) {

			String genotype;
			if (cnvs.containsKey(sample)) {
				if (cnvs.get(sample).getCopyNumber() == 1 || cnvs.get(sample).getCopyNumber() == 3) {
					genotype = "0/1:" + df.format(cnvs.get(sample).getConfidence());
					ac++;
				} else {
					genotype = "1/1:" + df.format(cnvs.get(sample).getConfidence());
					ac+=2;
				}
			} else {
				genotype = "0/0:.";
			}
		
			genotypeList.add(genotype);
			
		}
		
		writeRecord(genotypeList, chr, start, end, ct, ac);
				
	}
	
	private void writeRecord(List<String> genotypes, String chr, int start, int stop, CopyType ct, int ac) throws IOException, OverlapError {

		String alt = new String(fastaRef.getSubsequenceAt(chr, (long) start, (long) start).getBases());
		Set<String> genes = rawCNVAnnotator.parseGenes(chr, start, stop);
		Set<String> pathogenic = rawCNVAnnotator.parsePathogenic(chr, start, stop, genes, ct);
				
		
		String geneString = Combine.combineSet(genes, "|");
		String pathString = deconvolutePathogenic(pathogenic); //We should realistically only have one term here!
		int len = stop - start;
		if (ct == CopyType.DEL) {
			vcfWriter.write(chr + "\t" + start + "\t.\t" + alt + "\t<CN0>\t.\t.\tEND=" + stop + ";SVLEN=" + len + ";AC=" + ac + ";GENE=" + geneString + ";PATH=" + pathString + "\tGT:WC\t" + Combine.combineList(genotypes, "\t") + "\n");
		} else {
			vcfWriter.write(chr + "\t" + start + "\t.\t" + alt + "\t<DUP>\t.\t.\tEND=" + stop + ";SVLEN=" + len + ";AC=" + ac + ";GENE=" + geneString + ";PATH=" + pathString + "\tGT:WC\t" + Combine.combineList(genotypes, "\t") + "\n");
		}
		vcfWriter.flush();
			
	}

	private String deconvolutePathogenic(Set<String> pathogenic) {
		
		//This just sets up a heirarchy for commonly overlapped calls
		if (pathogenic.size() == 1) {
			return(pathogenic.iterator().next());
		} else {
			if (pathogenic.contains("LargeDel")) {
				return("LargeDel");
			} else if (pathogenic.contains("LargeDup")) {
				return("LargeDup");
			} else if (pathogenic.contains("15q13.3del")) {
				return("15q13.3del");
			} else if (pathogenic.contains("15q13.3dup")) {
				return("15q13.3dup");
			} else if (pathogenic.contains("2q13del")) {
				return("2q13del");
			} else if (pathogenic.contains("2q13dup")) {
				return("2q13dup");
			} else {
				return(null);
			}
		}
		
		
	}
	
}
