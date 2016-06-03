package projects.gemoma;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.EcNumber;

import de.jstacs.DataType;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.JstacsTool.ResultEntry;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class enables to extract introns from SAM files, which might be used to define splice sites in GeMoMa.
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
		
		public static void addIntrons(String[] samLine, Stranded stranded, List<Intron> introns){
			
			int start = Integer.parseInt(samLine[3]);
			
			String cigar = samLine[5];
			if(!cigar.contains("N")){
				return;
			}
			
			int bitflag = Integer.parseInt(samLine[1]);
			
			startOffs.clear();
			lens.clear();
			getOffset(cigar,startOffs,lens);
			Strand strand = getStrand(bitflag, stranded);
			
			for(int i=0;i<startOffs.length();i++){
				Intron in = new Intron( start+startOffs.get(i), lens.get(i),strand );
				introns.add(in);
			}
			
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
		/*TODO uncomment
		Stranded stranded = (Stranded) parameters.getParameterAt(0).getValue();
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterAt(1).getValue();
				
		int i=0;
		String str = null;
		
		BufferedReader reader;
		for( int k = 0; k < eps.getNumberOfParameters(); k++ ) {
			String fName = ((ParameterSet)eps.getParameterAt(k).getValue()).getParameterAt(0).getValue().toString();
			//System.out.println(fName);

			reader = new BufferedReader(new FileReader(fName));
			while( (str = reader.readLine()) != null ){
				if( str.charAt(0) == '@' ) continue;
				
				String[] parts = str.split("\t");
				String chrom = parts[2];
				
				ArrayList<Intron> introns = intronMap.get(chrom);
				if( introns == null ) {
					introns = new ArrayList<ExtractIntrons.Intron>();
					intronMap.put(chrom, introns);
				}
				
				Intron.addIntrons(parts, stranded, introns);
				
				if(i % 1000000 == 0){
					System.err.println(i);
				}
				i++;
				
			}
			reader.close();
		}
		/**/
		File out = File.createTempFile("intron_gff", "_GeMoMa.temp", new File("."));
		out.deleteOnExit(); 

		Iterator<String> it = intronMap.keySet().iterator();
		
		SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(new FileOutputStream(out));
		while(it.hasNext()){
			String chrom = it.next();
			List<Intron> introns = intronMap.get(chrom);
			introns = count(introns);
			print(chrom,introns,sos);
		}		
		sos.close();

		return new ToolResult("", "", null, new ResultSet(new TextResult("predicted annotation", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff", getToolName(), null, true)), parameters, getToolName(), new Date());
	}

	private static void print(String chrom, List<Intron> introns,SafeOutputStream sos) throws IOException{
		Iterator<Intron> it = introns.iterator();
		while(it.hasNext()){
			Intron in = it.next();
			sos.writeln(chrom+"\tsam\tintron\t"+in.getStart()+"\t"+in.getEnd()+"\t"+in.getCount()+"\t"+in.getStrand()+"\t.\t.");
		}
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
				/*
				new SimpleParameterSet(
					new EnumParameter(Stranded.class, "Defines whether the reads are stranded", true),
					new ParameterSetContainer( new ExpandableParameterSet( new SimpleParameterSet(		
							new FileParameter( "mapped reads file", "a SAM file containing the mapped reads", "sam",  true )
						), "mapped reads", "", 1 ) )
				);/**/
					
				//complex
				new ExpandableParameterSet( new SimpleParameterSet(
					new EnumParameter(Stranded.class, "Defines whether the reads are stranded", true),
					new FileParameter( "mapped reads file", "a SAM file containing the mapped reads", "sam",  true )
				), "mapped reads", "", 1 );			
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
		return "extracts introns from SAM files";
	}
	
	public String getHelpText() {
		return 
			"**What it does**\n\nThis tools extracts introns from SAM files obtained from RNA-seq analysis, that can be used to define splice sites in GeMoMa.\n\n"
				+ "**References**\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa or contact jens.keilwagen@julius-kuehn.de.\n";

	}
	
	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "introns"),
		};
	}
}
