package projects.gemoma;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * This class enables to extract introns from BAM/SAM files, which might be used to define splice sites in GeMoMa.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see GeMoMa
 */
public class ExtractIntrons implements JstacsTool {

	private enum Stranded{
		NO,
		FIRSTFWD,
		SECONDFWD
	}
	
	private enum Strand{
		FWD,
		REV,
		UNK;
		
		public String toString(){
			switch(this){
			case FWD: return "+";
			case REV: return "-";
			default: return ".";
			}
		}
	}
	
	private static class Intron implements Comparable<Intron>{
		
		private static IntList startOffs = new IntList();
		private static IntList lens = new IntList();
		
		private int start;
		private int end;
		private Strand strand;
		private int count;
		
		public int getCount() {
			return count;
		}

		public void setCount(int count) {
			this.count = count;
		}

		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}

		public Strand getStrand() {
			return strand;
		}

		private Intron(int start, int len, Strand strand){			
			this.start = start;
			this.end = start+len;
			this.strand = strand;
			this.count = 0;
		}
		
		public static boolean addIntrons(SAMRecord record, Stranded stranded, List<Intron> introns){
			
			int start = record.getAlignmentStart();
			
			String cigar = record.getCigarString();
			if(!cigar.contains("N")){
				return false;
			}
			
			int bitflag = record.getFlags();
			
			startOffs.clear();
			lens.clear();
			getOffset(cigar,startOffs,lens);
			Strand strand = getStrand(bitflag, stranded);
			
			for(int i=0;i<startOffs.length();i++){
				Intron in = new Intron( start+startOffs.get(i), lens.get(i),strand );
				introns.add(in);
			}
			return startOffs.length()>0;
		}
	
		private static Strand getStrand(int bitflag, Stranded stranded) {
			if(stranded == Stranded.NO){
				return Strand.UNK;
			}else if(stranded == Stranded.FIRSTFWD){
				
				if( (bitflag & 128) == 128 ){// second in pair
					if((bitflag & 16) == 16 ){
						return Strand.FWD;
					}else{
						return Strand.REV;
					}
				}else{// first in pair or unpaired
					if((bitflag & 16) == 16 ){
						return Strand.REV;
					}else{
						return Strand.FWD;
					}
				}
				
			}else{
				
				if( (bitflag & 128) == 128 ){// second in pair
					if((bitflag & 16) == 16 ){
						return Strand.REV;
					}else{
						return Strand.FWD;
					}
				}else{// first in pair or unpaired
					if((bitflag & 16) == 16 ){
						return Strand.FWD;
					}else{
						return Strand.REV;
					}
				}
			}
		}

		private static Pattern p = Pattern.compile( "[0-9]+(M|D|N)" );
		
		private static void getOffset(String cigar, IntList startOffs, IntList lens){
			Matcher m = p.matcher(cigar);
			
			int off = 0;
			while( m.find() ){
				int len = Integer.parseInt( cigar.substring( m.start(),m.end()-1 ) );
				if(m.group().endsWith("N")){
					startOffs.add(off);
					lens.add(len);
				}
				off += len;
			}
		}

		@Override
		public int compareTo(Intron o) {
			int comp = this.strand.compareTo(o.strand);
			if(comp == 0){
				comp = Integer.compare(this.start, o.start);
			}
			if(comp == 0){
				comp = Integer.compare(this.end, o.end);
			}
			return comp;
		}
		
	}
	
	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		HashMap<String, ArrayList<Intron>> intronMap = new HashMap<String, ArrayList<Intron>>();

		Stranded stranded = (Stranded) parameters.getParameterAt(0).getValue();
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterAt(1).getValue();
		ValidationStringency stringency = (ValidationStringency) parameters.getParameterAt(2).getValue();
		
		long i=0, s=0;
		
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( stringency );//important for unmapped reads
		
		//BufferedReader reader;
		int corrupt = 0;
		for( int k = 0; k < eps.getNumberOfParameters(); k++ ) {
			String fName = ((ParameterSet)eps.getParameterAt(k).getValue()).getParameterAt(0).getValue().toString();
			protocol.append(fName+" " + (new Date()) + "\n");
			
			long a = 0, b = 0;
			SamReader sr = srf.open(new File(fName));
			SAMRecordIterator samIt = sr.iterator();
			try {
				while(samIt.hasNext()){
					SAMRecord rec = samIt.next();
					String chrom = rec.getReferenceName();	
					ArrayList<Intron> introns = intronMap.get(chrom);
					if( introns == null ) {
						introns = new ArrayList<ExtractIntrons.Intron>();
						intronMap.put(chrom, introns);
					}
					
					if( Intron.addIntrons(rec, stranded, introns) ) {
						s++;
						b++;
					}
					
					if(i % 1000000 == 0){
						protocol.append(i+"\n");
					}
					i++;
					a++;
					
				}
			} catch( Exception e ) {
				//even if the file is broken take all information before it breaks and then write a message
				e.printStackTrace();
				corrupt++;
			} finally {
				samIt.close();
				sr.close();
				
				protocol.append("statistics for " + fName +"\n");
				protocol.append("#reads: " + a +"\n");
				protocol.append("#split reads: " + b +"\n");
			}
		}
		
		File out = File.createTempFile("intron_gff", "_GeMoMa.temp", new File("."));
		out.deleteOnExit(); 

		Iterator<String> it = intronMap.keySet().iterator();
		
		SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(new FileOutputStream(out));
		int intronNum = 0;
		while(it.hasNext()){
			String chrom = it.next();
			List<Intron> introns = intronMap.get(chrom);
			introns = count(introns);
			intronNum += print(chrom,introns,sos);
		}		
		sos.close();

		protocol.append("\noverall statistics\n");
		protocol.append("#corrupt files: " + corrupt +"\n");
		protocol.append("#reads: " + i +"\n");
		protocol.append("#split reads: " + s +"\n");
		protocol.append("#introns: " + intronNum +"\n");
		
		return new ToolResult("", "", null, new ResultSet(new TextResult("introns", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff", getToolName(), null, true)), parameters, getToolName(), new Date());
	}

	private static int print(String chrom, List<Intron> introns,SafeOutputStream sos) throws IOException{
		Iterator<Intron> it = introns.iterator();
		int i = 0;
		while(it.hasNext()){
			Intron in = it.next();
			sos.writeln(chrom+"\tsam\tintron\t"+in.getStart()+"\t"+in.getEnd()+"\t"+in.getCount()+"\t"+in.getStrand()+"\t.\t.");
			i++;
		}
		return i;
	}
	
	private static List<Intron> count(List<Intron> introns) {
		
		Collections.sort(introns);
		
		ArrayList<Intron> agg = new ArrayList<ExtractIntrons.Intron>();
		
		if(introns.size() == 0){
			return agg;
		}
		
		Intron last = introns.get(0);
		int n=1;
		for(int i=1;i<introns.size();i++){
			Intron curr = introns.get(i);
			if(last.compareTo(curr) == 0){
				n++;
			}else{
				last.setCount(n);
				agg.add(last);
				last = curr;
				n = 1;
			}
		}
		
		last.setCount(n);
		agg.add(last);
		
		return agg;
	}

	@Override
	public ParameterSet getToolParameters() {
		try{
			return
				//TODO for testing: which is better for the end-users 
				//simple
				
				new SimpleParameterSet(
					new EnumParameter(Stranded.class, "Defines whether the reads are stranded", true),
					new ParameterSetContainer( new ExpandableParameterSet( new SimpleParameterSet(		
							new FileParameter( "mapped reads file", "BAM/SAM files containing the mapped reads", "bam,sam",  true )
						), "mapped reads", "", 1 ) ),
					new EnumParameter(ValidationStringency.class, "Defines how strict to be when reading a SAM or BAM, beyond bare minimum validation.", true, ValidationStringency.LENIENT.name() )
				);
				/*					
				//complex
				new ExpandableParameterSet( new SimpleParameterSet(
					new EnumParameter(Stranded.class, "Defines whether the reads are stranded", true),
					new FileParameter( "mapped reads file", "a SAM file containing the mapped reads", "sam",  true )
				), "mapped reads", "", 1 );/**/			
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return getShortName();
	}

	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public String getShortName() {
		return "ExtractIntrons";
	}

	@Override
	public String getDescription() {
		return "extracts introns from BAM/SAM files";
	}
	
	public String getHelpText() {
		return 
			"**What it does**\n\nThis tools extracts introns from BAM/SAM files obtained from RNA-seq analysis, that can be used to define splice sites in GeMoMa.\n\n"
				+ "**References**\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa or contact jens.keilwagen@julius-kuehn.de.\n";

	}
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "introns"),
		};
	}
}
