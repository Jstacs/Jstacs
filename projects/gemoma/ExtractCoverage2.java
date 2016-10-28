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

import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.SafeOutputStream;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import projects.gemoma.ExtractIntrons.Stranded;

public class ExtractCoverage2 implements JstacsTool {

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
	
	
	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new ExtractCoverage2());
		
		cli.run(args);
	}
	
	
	@Override
	public ParameterSet getToolParameters() {
		
		try {
			return new SimpleParameterSet(
						new EnumParameter(Stranded.class, "Defines whether the reads are stranded", true),
						new ParameterSetContainer( new ExpandableParameterSet( new SimpleParameterSet(		
								new FileParameter( "mapped reads file", "BAM/SAM files containing the mapped reads", "bam,sam",  true )
							), "mapped reads", "", 1 ) )
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
				 
		
		
		//int i=0;
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		
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
		
		
		
		File outFwd = File.createTempFile("coveragefwd_bedgraph", "_GeMoMa.temp", new File("."));
		outFwd.deleteOnExit(); 
		
		SafeOutputStream sosFwd = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outFwd));
		sosFwd.writeln("track type=bedgraph");
		
		File outRev = File.createTempFile("coveragerev_bedgraph", "_GeMoMa.temp", new File("."));
		outRev.deleteOnExit(); 
		
		SafeOutputStream sosRev = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outRev));
		sosRev.writeln("track type=bedgraph");
		
		
		
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
		while(true){
			
			int wm = whichMin(curr,chr);
			
			if(wm > -1){
			
				SAMRecord rec = curr[wm];
				if( !rec.isSecondaryOrSupplementary() ) {//TODO?
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
				if(its[wm].hasNext()){
					curr[wm] = its[wm].next();
				}else{
					curr[wm] = null;
				}
			
			}else{
				
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
				
				chr = null;
				for(int k=0;k<curr.length;k++){
					if(curr[k] != null){
						firstChrs[k] = curr[k].getReferenceName();
					}else{
						firstChrs[k] = null;
					}
				}
				Arrays.sort(firstChrs,scomp);
				chr = firstChrs[0];
				mapFwd.clear();
				mapRev.clear();
				currPos = 0;
				if(chr == null){
					break;
				}
			}
		}

		sosFwd.close();
		sosRev.close();
		return new ToolResult("", "", null, new ResultSet(new Result[]{new TextResult("coverage forward", "Result", new FileParameter.FileRepresentation(outFwd.getAbsolutePath()), "bedgraph", getToolName(), null, true),
				new TextResult("coverage reverse", "Result", new FileParameter.FileRepresentation(outRev.getAbsolutePath()), "bedgraph", getToolName(), null, true)}), parameters, getToolName(), new Date());
		
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
		return "Extract Coverage";
	}

	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public String getShortName() {
		return "coverage";
	}

	@Override
	public String getDescription() {
		return "extract coverage from SAM/BAM";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	
		
	
	
}
