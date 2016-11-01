package projects.gemoma;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.JstacsTool.ResultEntry;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFileHeader.SortOrder;

/**
 * This class enables to extract the coverage per strand and introns from BAM/SAM files, which might be used in GeMoMa.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ExtractRNAseqEvidence implements JstacsTool {

	public static class BedgraphEntry{
		
		private String chr;
		private int start;
		private int end;
		private int value;
		
		public BedgraphEntry(String chr, int start, int end, int value){
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.value = value;
		}
		
	}
	
	public enum Stranded{
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
	
	
	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new ExtractRNAseqEvidence());
		
		cli.run(args);
	}
	
	
	@Override
	public ParameterSet getToolParameters() {
		
		try {
			return new SimpleParameterSet(
						new EnumParameter(Stranded.class, "Defines whether the reads are stranded", true),
						new ParameterSetContainer( new ExpandableParameterSet( new SimpleParameterSet(		
								new FileParameter( "mapped reads file", "BAM/SAM files containing the mapped reads", "bam,sam",  true )
							), "mapped reads", "", 1 ) ),
						new EnumParameter(ValidationStringency.class, "Defines how strict to be when reading a SAM or BAM, beyond bare minimum validation.", true, ValidationStringency.LENIENT.name() ),
						new SimpleParameter(DataType.BOOLEAN,"use secondary alignments", "allows to filter flags in the SAM or BAM", true, true),
						new SimpleParameter(DataType.BOOLEAN,"coverage output", "allows to output the coverage", true, false)
					);
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException();
		}
		
	}

	@Override
	public ToolResult run(ParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		Stranded stranded = (Stranded) parameters.getParameterAt(0).getValue();
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterAt(1).getValue();
		ValidationStringency stringency = (ValidationStringency) parameters.getParameterAt(2).getValue();
		boolean sa = (Boolean) parameters.getParameterAt(3).getValue();
		boolean coverage = (Boolean) parameters.getParameterAt(4).getValue();
		
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( stringency );//important for unmapped reads
		
		SAMRecordIterator[] its = new SAMRecordIterator[eps.getNumberOfParameters()];
		SAMRecord[] curr = new SAMRecord[eps.getNumberOfParameters()];
		String[] firstChrs = new String[eps.getNumberOfParameters()];
		for( int k = 0; k < eps.getNumberOfParameters(); k++ ) {
			String fName = ((ParameterSet)eps.getParameterAt(k).getValue()).getParameterAt(0).getValue().toString();
			SamReader sr = srf.open(new File(fName));
			
			SAMRecordIterator samIt = sr.iterator();
			samIt = samIt.assertSorted(SortOrder.coordinate);
			its[k] = samIt;
			if(its[k].hasNext()){
				curr[k] = its[k].next();
				firstChrs[k] = curr[k].getReferenceName();
			}
		}
		
		
		SafeOutputStream sosFwd, sosRev;
		File outFwd = null, outRev = null;
		if( coverage ) {
			outFwd = File.createTempFile("coveragefwd_bedgraph", "_GeMoMa.temp", new File("."));
			outFwd.deleteOnExit(); 
		
			sosFwd = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outFwd));
			
			outRev = File.createTempFile("coveragerev_bedgraph", "_GeMoMa.temp", new File("."));
			outRev.deleteOnExit();
			
			sosRev = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outRev));
		} else {
			sosFwd = SafeOutputStream.getSafeOutputStream(null);
			sosRev = SafeOutputStream.getSafeOutputStream(null);
		}

		sosFwd.writeln("track type=bedgraph");
		sosRev.writeln("track type=bedgraph");
		
		File outInt = File.createTempFile("intron_gff", "_GeMoMa.temp", new File("."));
		outInt.deleteOnExit(); 
		SafeOutputStream sosInt = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outInt));	
		
		
		Comparator<String> scomp = new Comparator<String>() {
			
			@Override
			public int compare(String o1, String o2) {
				if(o1 == null && o2 == null){
					return 0;
				}else if(o1 == null){
					return 1;
				}else if(o2 == null){
					return -1;
				}else{
					return o1.compareTo(o2);
				}
			}
		};
		
		Arrays.sort(firstChrs,scomp);
		
		HashMap<Integer, int[]> mapFwd = new HashMap<Integer, int[]>();
		HashMap<Integer, int[]> mapRev = new HashMap<Integer, int[]>();
		
		BedgraphEntry previousFwd = null;
		BedgraphEntry previousRev = null;
		
		String chr = firstChrs[0];
		int currPos = 0;
		ArrayList<Intron> introns = new ArrayList<Intron>();
		long splits = 0, intronNum = 0;
		long[] split = new long[its.length];
		boolean[] corrupt = new boolean[its.length];
		Arrays.fill(corrupt, false);
		long i = 0;
		while(true){
			
			int wm = whichMin(curr,chr);
			
			if(wm > -1){
			
				SAMRecord rec = curr[wm];
				if( sa || !rec.isSecondaryOrSupplementary() ) {
					if( coverage ) {
						boolean isNeg = rec.getReadNegativeStrandFlag();
						boolean isFirst = !rec.getReadPairedFlag() || rec.getFirstOfPairFlag();
						boolean countAsFwd = true;
						if(stranded == Stranded.FIRSTFWD){
							countAsFwd = (isFirst && !isNeg)||(!isFirst && isNeg);
						}else if(stranded == Stranded.SECONDFWD){
							countAsFwd = (isFirst && isNeg)||(!isFirst && !isNeg);
						}
						
						
						
						int recStart = rec.getStart();
						while(currPos < recStart){
							if(mapFwd.containsKey(currPos)){
								previousFwd = write(previousFwd,chr,mapFwd,currPos,sosFwd);
								mapFwd.remove(currPos);
							}
							if(mapRev.containsKey(currPos)){
								previousRev = write(previousRev,chr,mapRev,currPos,sosRev);
								mapRev.remove(currPos);
							}
							currPos++;
						}
						
						HashMap<Integer, int[]> map = countAsFwd ? mapFwd : mapRev;
						List<AlignmentBlock> blocks = rec.getAlignmentBlocks();
						Iterator<AlignmentBlock> blockIt = blocks.iterator();
						while(blockIt.hasNext()){
							AlignmentBlock block = blockIt.next();
							int start = block.getReferenceStart();
							int len = block.getLength();
							for(int k=0;k<len;k++){
								if(!map.containsKey(start+k)){
									map.put(start+k, new int[1]);
								}
								map.get(start+k)[0]++;
							}
						}
					}
					if( Intron.addIntrons(rec, stranded, introns) ) {
						splits++;
						split[wm]++;
					}
				}
				if(its[wm].hasNext()){
					try {
						curr[wm] = its[wm].next();
					} catch( Exception e ) {
						//even if the file is broken take all information before it breaks and then write a message
						e.printStackTrace();
						corrupt[wm]=true;
						curr[wm]=null;
						protocol.append("corrupt file: " + ((ParameterSet)eps.getParameterAt(wm).getValue()).getParameterAt(0).getValue().toString() + "\n" );
					}
				}else{
					curr[wm] = null;
				}
			
			}else{
				if( coverage ) {
					int temp = currPos;
					while(mapFwd.size()>0){
						if(mapFwd.containsKey(temp)){
							previousFwd = write(previousFwd,chr,mapFwd,temp,sosFwd);
							mapFwd.remove(temp);
						}
						temp++;
					}
					temp = currPos;
					while(mapRev.size()>0){
						if(mapRev.containsKey(temp)){
							previousRev = write(previousRev,chr,mapRev,temp,sosRev);
							mapRev.remove(temp);
						}
						temp++;
					}
				}
				for(int k=0;k<curr.length;k++){
					if(curr[k] != null){
						firstChrs[k] = curr[k].getReferenceName();
					}else{
						firstChrs[k] = null;
					}
				}
				Arrays.sort(firstChrs,scomp);
				
				if( firstChrs[0] == null || !chr.equals(firstChrs[0]) ) {
					introns = count(introns);
					intronNum += print(chr,introns,sosInt);
					introns.clear();
				}
				
				chr = firstChrs[0];
				mapFwd.clear();
				mapRev.clear();
				currPos = 0;
				if(firstChrs[0] == null){
					break;
				}
				
				
			}
			
			i++;
			if(i % 1000000 == 0){
				protocol.append(i+"\n");
			}
		}
		sosInt.close();
		sosFwd.close();
		sosRev.close();
		
		protocol.append("\nfile statistics:\n");
		int c = 0;
		for( int k = 0; k < its.length; k++ ) {
			its[k].close();
			protocol.append( split[k] + "\t" + corrupt[k] + "\t" + ((ParameterSet)eps.getParameterAt(k).getValue()).getParameterAt(0).getValue().toString() + "\n" );
			c += (corrupt[k]?1:0);
		}
		protocol.append("\noverall statistics:\n");
		protocol.append("#files:\t" + its.length + "\n");
		protocol.append("#corrupt files:\t" + c + "\n");
		protocol.append("#reads:\t" + i + "\n");
		protocol.append("#split reads:\t" + splits + "\n");
		protocol.append("#introns:\t" + intronNum + "\n");
		
		Result[] res = new Result[coverage?3:1];
		res[0] = new TextResult("introns", "Result", new FileParameter.FileRepresentation(outInt.getAbsolutePath()), "gff", getToolName(), null, true);
		if( coverage ) {
			res[1] = new TextResult("coverage forward", "Result", new FileParameter.FileRepresentation(outFwd.getAbsolutePath()), "bedgraph", getToolName(), null, true);
			res[2] = new TextResult("coverage reverse", "Result", new FileParameter.FileRepresentation(outRev.getAbsolutePath()), "bedgraph", getToolName(), null, true);
		}
		
		return new ToolResult("", "", null, new ResultSet( res ), parameters, getToolName(), new Date());
		
	}
	
	private static long print(String chrom, ArrayList<Intron> introns,SafeOutputStream sos) throws IOException{
		Iterator<Intron> it = introns.iterator();
		long i = 0;
		while(it.hasNext()){
			Intron in = it.next();
			sos.writeln(chrom+"\tsam\tintron\t"+in.getStart()+"\t"+in.getEnd()+"\t"+in.getCount()+"\t"+in.getStrand()+"\t.\t.");
			i++;
		}
		return i;
	}
	
	private static ArrayList<Intron> count(ArrayList<Intron> introns) {
		Collections.sort(introns);
		ArrayList<Intron> agg = new ArrayList<Intron>();
		
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
	

	private BedgraphEntry write(BedgraphEntry previous, String chr, HashMap<Integer, int[]> map, int currPos, SafeOutputStream safeOutputStream) throws IOException {
		int[] temp = map.get(currPos);
		
		
		if(temp != null){
			if( previous != null && currPos == previous.end && temp[0] == previous.value && chr.equals(previous.chr) ){
				previous.end = currPos+1;
				return previous;
			}else{
				if(previous != null){
					safeOutputStream.writeln(previous.chr+"\t"+previous.start+"\t"+(previous.end)+"\t"+previous.value);
				}
				return new BedgraphEntry(chr, currPos, currPos+1, temp[0]);
			}
		}else if(previous != null){
			safeOutputStream.writeln(previous.chr+"\t"+previous.start+"\t"+(previous.end)+"\t"+previous.value);
		}
		return null;
	}

	private int whichMin(SAMRecord[] curr, String chr) {
		int res = -1;
		for(int i=0;i<curr.length;i++){
			if(curr[i] != null && chr.equals(curr[i].getReferenceName())){
				if(res == -1 || curr[i].getStart() < curr[res].getStart()){
					res = i;
				}
			}
		}
		return res;
	}

	@Override
	public String getToolName() {
		return "Extract RNA-seq Evidence";
	}

	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public String getShortName() {
		return "ERE";
	}

	@Override
	public String getDescription() {
		return "extract introns and coverage from SAM/BAM";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "introns")
		};
	}
}
