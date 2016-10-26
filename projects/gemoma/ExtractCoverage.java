package projects.gemoma;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collections;
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
import projects.gemoma.ExtractIntrons.Stranded;

public class ExtractCoverage implements JstacsTool {

	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new ExtractCoverage());
		
		cli.run(args);
	}
	
	
	public static class Read implements Comparable<Read>{
		
		private int[] starts;
		private int[] lengths;
		private boolean countAsFwd;
		
		public Read(List<AlignmentBlock> blocks, boolean isNeg, boolean isFirst, Stranded stranded){
			starts = new int[blocks.size()];
			lengths = new int[blocks.size()];
			Iterator<AlignmentBlock> it = blocks.iterator();
			int i=0;
			while(it.hasNext()){
				AlignmentBlock ab = it.next();
				starts[i] = ab.getReferenceStart();
				lengths[i] = ab.getLength();
				i++;
			}
			
			if(stranded == Stranded.FIRSTFWD){
				countAsFwd = (isFirst && !isNeg)||(!isFirst && isNeg);
			}else if(stranded == Stranded.SECONDFWD){
				countAsFwd = (isFirst && isNeg)||(!isFirst && !isNeg);
			}else{
				countAsFwd = true;
			}
		}

		@Override
		public int compareTo(Read o) {
			return Integer.compare(starts[0], o.starts[0]);
		}

		public int getStart() {
			return starts[0];
		}

		public int getEnd() {
			return starts[starts.length-1] + lengths[lengths.length-1];
		}

		public boolean overlaps(int currPos) {
			int idx = 0;
			while( idx < starts.length && starts[idx]+lengths[idx]<=currPos ){
				idx++;
			}
			return idx < starts.length && starts[idx] <= currPos;
		}
		
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
				
		
		HashMap<String,ArrayList<Read>> map = new HashMap<String,ArrayList<Read>>(); 
		
		
		int i=0;
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		
		for( int k = 0; k < eps.getNumberOfParameters(); k++ ) {
			String fName = ((ParameterSet)eps.getParameterAt(k).getValue()).getParameterAt(0).getValue().toString();
					SamReader sr = srf.open(new File(fName));
			
			SAMRecordIterator samIt = sr.iterator();		
			while(samIt.hasNext()){
				SAMRecord rec = samIt.next();
				String chr = rec.getReferenceName();
				List<AlignmentBlock> blocks = rec.getAlignmentBlocks();
				boolean isNeg = rec.getReadNegativeStrandFlag();
				boolean isFirst = !rec.getReadPairedFlag() || rec.getFirstOfPairFlag();
				Read r = new Read(blocks,isNeg,isFirst,stranded);
				
				if(!map.containsKey(chr)){
					map.put(chr, new ArrayList<Read>());
				}
				map.get(chr).add(r);	
			}
		}
		
		Iterator<String> keyIt = map.keySet().iterator();
		
		File outFwd = File.createTempFile("coveragefwd_bedgraph", "_GeMoMa.temp", new File("."));
		outFwd.deleteOnExit(); 
		
		SafeOutputStream sosFwd = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outFwd));
		sosFwd.writeln("track type=bedgraph");
		
		File outRev = File.createTempFile("coveragerev_bedgraph", "_GeMoMa.temp", new File("."));
		outRev.deleteOnExit(); 
		
		SafeOutputStream sosRev = SafeOutputStream.getSafeOutputStream(new FileOutputStream(outRev));
		sosRev.writeln("track type=bedgraph");
		
		while(keyIt.hasNext()){
			String chr = keyIt.next();
			
			ArrayList<Read> list = map.get(chr);
			
			Collections.sort(list);
			
			int readOff = 0;
			int currPos = 0;
			
			while(readOff < list.size()){
				
				while(readOff<list.size() && list.get(readOff).getEnd()<=currPos){
					readOff++;
				}
				if(readOff == list.size()){
					break;
				}
				while(list.get(readOff).getStart()>currPos){
					currPos++;
				}
				
				int temp = readOff;
				int countFwd = 0;
				int countRev = 0;
				while(temp < list.size() && list.get(temp).getStart()<=currPos){
					Read r = list.get(temp);
					if(r.countAsFwd){
						if(r.overlaps(currPos)){
							countFwd++;
						}
					}else{
						if(r.overlaps(currPos)){
							countRev++;
						}
					}
					temp++;
				}
				if(countFwd > 0){
					sosFwd.writeln(chr+"\t"+currPos+"\t"+(currPos+1)+"\t"+countFwd);
				}
				if(countRev > 0){
					sosRev.writeln(chr+"\t"+currPos+"\t"+(currPos+1)+"\t"+countRev);
				}
				currPos++;
				
			}
			

			
		}

		sosFwd.close();
		sosRev.close();
		return new ToolResult("", "", null, new ResultSet(new Result[]{new TextResult("coverage forward", "Result", new FileParameter.FileRepresentation(outFwd.getAbsolutePath()), "bedgraph", getToolName(), null, true),
				new TextResult("coverage reverse", "Result", new FileParameter.FileRepresentation(outRev.getAbsolutePath()), "bedgraph", getToolName(), null, true)}), parameters, getToolName(), new Date());
		
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
