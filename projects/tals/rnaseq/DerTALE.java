package projects.tals.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
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
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class DerTALE implements JstacsTool {

	private enum Compare{
		EXTREMES,
		MEDIAN,
		MEAN
	}
	
	private static class Entry{
		private String chr;
		private int pos;
		private String strand;
		
		private LinkedList<double[]> positives;
		private LinkedList<double[]> negatives;
		
		
		public static String getHeader(){
			return "Chr\tPosition\tStrand";
		}
		
		public String toString(){
			return chr+"\t"+pos+"\t"+strand;
		}
		


		public void addPositives(double[] positives) {
			if(this.positives == null){
				this.positives = new LinkedList<>();
			}else{
				if(this.positives.getFirst().length != positives.length){
					throw new RuntimeException();
				}
			}
			this.positives.add(positives);
		}
		
		public void addNegatives(double[] negatives) {
			if(this.negatives == null){
				this.negatives = new LinkedList<>();
			}else{
				if(this.negatives.getFirst().length != negatives.length){
					throw new RuntimeException();
				}
			}
			this.negatives.add(negatives);
		}

		
		public Entry(String chr, int pos, String strand){
			this.chr = chr;
			this.pos = pos;
			this.strand = strand;
		}

		public String getChr() {
			return chr;
		}

		public int getPos() {
			return pos;
		}

		public String getStrand() {
			return strand;
		}

		
		private static double[] getSmoothedProfile(double[] profile, int windowWidth){
			double currPos = ToolBox.sum(0, windowWidth, profile);
			
			double[] res = new double[profile.length-windowWidth];
			for(int i=windowWidth;i<profile.length;i++){
				res[i-windowWidth] = currPos;
				currPos += (profile[i] - profile[i-windowWidth]);
				//currPos = ToolBox.median(i-windowWidth+1,i+1,profile);
			}
			return res;
		}
		
		private static double[][] getSmoothedProfiles(LinkedList<double[]> profs, int windowWidth){
			double[][] dProfs = new double[profs.size()][];
			for(int i=0;i<dProfs.length;i++){
				dProfs[i] = getSmoothedProfile(profs.get(i),windowWidth);
			}
			return dProfs;
		}
		

		public Pair<double[],String> getProfileResult(String[] header,int windowWidth, double t, Compare comp) {
			double[][] posProfs = getSmoothedProfiles(positives, windowWidth);
			double[][] negProfs = getSmoothedProfiles(negatives, windowWidth);
			
			int mid = positives.getFirst().length/2; //3000
			
			int relIdx = -1;
			double max = Double.NEGATIVE_INFINITY;
			
			int sumAboveThreshold=0;
			int[] posAboveThreshold=new int[posProfs[0].length];
			Arrays.fill(posAboveThreshold, Integer.MIN_VALUE);
			
			
			
			StringBuffer sb = new StringBuffer();
			sb.append("Position");
			for(int i=0;i<header.length;i++){
				sb.append("\t"+header[i].replaceAll(".*/SRR", "SRR").replace("/", "_"));
				
			}
			sb.append("\tabove threshold\n");
			
			for(int i=0;i<windowWidth/2;i++){
			
				sb.append(this.pos-mid+i);
				for(int j=0;j<posProfs.length;j++){
					sb.append("\t"+positives.get(j)[i]);
				}
				for(int j=0;j<negProfs.length;j++){
					sb.append("\t"+negatives.get(j)[i]);
				}
				sb.append("\tNA\n");
			}
			
			
	
			for(int i=0;i<posProfs[0].length;i++){
				
				sb.append(this.pos-mid+i+windowWidth/2);
				for(int j=0;j<posProfs.length;j++){
					sb.append("\t"+positives.get(j)[i+windowWidth/2]);
				}
				for(int j=0;j<negProfs.length;j++){
					sb.append("\t"+negatives.get(j)[i+windowWidth/2]);
				}
				
				double minPos = Double.POSITIVE_INFINITY;
				double maxNeg = Double.NEGATIVE_INFINITY;
				if(comp==Compare.EXTREMES){
					for(int j=0;j<posProfs.length;j++){
						if(posProfs[j][i] < minPos){
							minPos = posProfs[j][i];
						}
					}

					for(int j=0;j<negProfs.length;j++){
						if(negProfs[j][i] > maxNeg){
							maxNeg = negProfs[j][i];
						}
					}
				}else{
					double[] temp = new double[posProfs.length];
					for(int j=0;j<posProfs.length;j++){
						temp[j] = posProfs[j][i];
						
					}
					if(comp == Compare.MEAN){
						minPos = ToolBox.mean(temp);
					}else{
						minPos = ToolBox.median(temp);
					}
					temp = new double[negProfs.length];
					for(int j=0;j<negProfs.length;j++){
						temp[j] = negProfs[j][i];
					}
					if(comp == Compare.MEAN){
						maxNeg = ToolBox.mean(temp);
					}else{
						maxNeg = ToolBox.median(temp);
					}
				}
				double rat = Math.log(minPos) - Math.log(maxNeg);
				
				if(rat > max){
					max = rat;
					relIdx = i+windowWidth/2-mid;
				}
				
				if(rat > t){   
					posAboveThreshold[i]=i+windowWidth/2-mid;
					sumAboveThreshold++;
					sb.append("\tTRUE\n");
				}else{
					sb.append("\tFALSE\n");
				}
			}
			boolean prevPovAbove=false;
			int maxStrechLength=0;
			int aktStrechAboveThreshold=0;
			ArrayList<Integer> strechLengths=new ArrayList<Integer>();
			int whichStrechSpans=-1;
			int countSpans=0;
			int spanWidth=50;
		
			for(int i=0;i<posProfs[0].length;i++){
				
				if((posAboveThreshold[i]<=spanWidth)&&(posAboveThreshold[i]>=-spanWidth)){
					countSpans++;
				}
				
				if(posAboveThreshold[i]>Integer.MIN_VALUE){
					prevPovAbove=true;
					aktStrechAboveThreshold++;
				}else {
					if(aktStrechAboveThreshold>maxStrechLength){
				
						maxStrechLength=aktStrechAboveThreshold;
						strechLengths.add(aktStrechAboveThreshold);
						if((whichStrechSpans==-1)&&(countSpans==(spanWidth*2+1))){
							whichStrechSpans=strechLengths.size()-1;
						}
						aktStrechAboveThreshold=0;
						prevPovAbove=false;
						
					}
					
				}
			}
			int minStrechLength=400;
			boolean secStrech=false;
			if(whichStrechSpans!=-1){
				for(int i=0;i<strechLengths.size();i++){
					if(i!=whichStrechSpans){
						if(strechLengths.get(i)>=minStrechLength){
							secStrech=true;
						}
					}
				}
			}
			
			
			for(int i=posProfs[0].length+windowWidth/2;i<posProfs[0].length+windowWidth;i++){
				sb.append(this.pos-mid+i);
				for(int j=0;j<posProfs.length;j++){
					sb.append("\t"+positives.get(j)[i]);
				}
				for(int j=0;j<negProfs.length;j++){
					sb.append("\t"+negatives.get(j)[i]);
				}
				sb.append("\tNA\n");
			}
			
			
			if((max > t)&&(maxStrechLength>=minStrechLength)){
				
				if((whichStrechSpans!=-1)&&(secStrech)){
					return new Pair<double[], String>(new double[]{max,relIdx}, sb.toString());
				}
				else if(whichStrechSpans==-1){
					
					return new Pair<double[], String>(new double[]{max,relIdx}, sb.toString());
				}else {
					return new Pair<double[], String>(new double[]{max,relIdx}, null);
				}
				
				
			}else{
				
				return new Pair<double[], String>(new double[]{max,relIdx}, null);
			}
			
			
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new DerTALE());
		
		cli.run(args);

	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		try {
			pars.add(new FileParameter("Predictions", "Predictions output file", "tsv,tabular", true));
			pars.add(new ParameterSetContainer( "Treatment","",new ExpandableParameterSet(new SimpleParameterSet(new FileParameter("Treatment BAM", "BAM file of mapped reads from treatment experiment. BAM file must have an index with additional extension .bai.", "bam", true)), "Treatment data", "")) );
			pars.add(new ParameterSetContainer( "Control","",new ExpandableParameterSet(new SimpleParameterSet(new FileParameter("Control BAM", "BAM file of mapped reads from control experiment. BAM file must have an index with additional extension .bai.", "bam", true)), "Control data", "")) );
		} catch (CloneNotSupportedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "Number of predictions", "Number of (top) predictions considered", true,100));
			pars.add(new SimpleParameter(DataType.INT, "Region width", "Number of bases around the predicted site", true,3000));
			pars.add(new SimpleParameter(DataType.INT, "Window width", "Width of the window considered for differential abundance", true,300));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Pseudo count","Pseudo count on the count profile",true,1.0));
			pars.add(new EnumParameter(Compare.class, "Measure for comparing replicates", true,"MEAN"));
			pars.add(new SimpleParameter(DataType.DOUBLE, "Threshold", "Threshold on the log differential abundance", true,1.0));
		} catch (ParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		return new ToolParameterSet(this.getShortName(), pars.toArray(new Parameter[0]));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		
		FileRepresentation predsFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents();
		
		
		String[] posFiles = getFiles((ExpandableParameterSet)parameters.getParameterAt(1).getValue());
		String[] negFiles = getFiles((ExpandableParameterSet)parameters.getParameterAt(2).getValue());
		String[] files = new String[posFiles.length+negFiles.length];
		System.arraycopy(posFiles, 0, files, 0, posFiles.length);
		System.arraycopy(negFiles, 0, files, posFiles.length, negFiles.length);
		
		
		int numTop = (int) parameters.getParameterAt(3).getValue(); //100
		int regionWidth = (int) parameters.getParameterAt(4).getValue(); //3000
		int windowWidth = (int) parameters.getParameterAt(5).getValue(); //300
		double pseudoCount = (double) parameters.getParameterAt(6).getValue(); //1.0
		Compare compare = (Compare)((EnumParameter)parameters.getParameterAt(7)).getValue(); //MEAN
		double t = (double) parameters.getParameterAt(8).getValue();
		
		LinkedList<Entry> predsList = new LinkedList<>(); //einlsen der Liste an Predictions für einen TAL
		BufferedReader br = new BufferedReader(new StringReader(predsFile.getContent())); //BXOR1_RVDs/Predicted_binding_sites_for_TalAD9.tsv
		//# Seq-ID	Position	Strand	Score	Sequence	Approx. p-value	RVDs	TALE
		String curr = null;
		if( (curr = br.readLine()) != null ){
			//targetFinder  //Sequence Name   Strand  Score   Start Position  Target Sequence Plus strand sequence
			if(curr.startsWith("options_used:")){
				while( (curr = br.readLine()) != null ){
					if(!((curr.startsWith("Best"))||(curr.startsWith("Sequence")))){
						curr=curr.replaceAll("\\s+", " ");
						String[] parts = curr.split(" ");
						String strand=parts[1];
			
						if(strand.equals("Plus")){
							strand="+";
						}else if (strand.equals("Minus")) {
							strand="-";
						}
						predsList.add(new Entry(parts[0], Integer.parseInt(parts[3]), strand));
						
					}
					if(predsList.size()>=numTop){
						break;
					}
				}
			}
			//Talvez //TODO Warum "-" in File? genomische Position???!!! //>TalAD4 Chr9    17.357  -4893450        1
			else if (curr.startsWith(">")) {
				curr=curr.replaceAll("\\s+", " ");
				String[] parts = curr.split(" ");
				String strand =parts[1];
				if(strand.contains("_revcom")){
					strand="-";
				}else{
					strand="+";
				}
				predsList.add(new Entry(parts[1].replace("_revcom", ""), Integer.parseInt(parts[3]), strand)); 
				
				while( (curr = br.readLine()) != null ){
					curr=curr.replaceAll("\\s+", " ");
					parts = curr.split(" ");
					strand =parts[1];
					if(strand.contains("_revcom")){
						strand="-";
					}else{
						strand="+";
					}
					predsList.add(new Entry(parts[1].replace("_revcom", ""), Integer.parseInt(parts[3]), strand));
					
					if(predsList.size()>=numTop){
						break;
					}
				}
				//PrediTAL
			}else if (curr.startsWith("#")) {
				while( (curr = br.readLine()) != null ){
					
					if(!curr.startsWith("#")){
						String[] parts = curr.split("\t");
						predsList.add(new Entry(parts[0], Integer.parseInt(parts[1]), parts[2])); 
					}
					if(predsList.size()>=numTop){
						break;
					}
				}
			}else {
				throw new RuntimeException("File format of the prediction file does not look like the output of PrediTALE, Talvez or Target Finder!");
			}
		}
	
		
		br.close();
		
		
		
		
		
		SamReaderFactory srf = SamReaderFactory.makeDefault();
		srf.validationStringency( ValidationStringency.SILENT );
		
		for( int k = 0; k < files.length; k++ ) {
			String fName = files[k];
			
			
			for(Entry en : predsList){
				
				if(!new File(fName+".bai").exists()){
					throw new RuntimeException("No index found for file "+fName+". The index must be in the same directory as the specified BAM file with filename "+fName+".bai.");
				}
				
				SamReader sr = srf.open(SamInputResource.of(new File(fName)).index(new File(fName+".bai")));
				
				


				int start = en.getPos()-regionWidth; //3000 positionen davor und 3000 Positionen danach
				int end = en.getPos()+regionWidth;
				SAMRecordIterator sri = sr.query(en.chr, start, end, true);
				
				BAMIndex index = sr.indexing().getIndex();
				double count = 0;
		        for (int i = 0; i < sr.getFileHeader().getSequenceDictionary().size(); i++) {
		            BAMIndexMetaData meta = index.getMetaData(i);
		            count += meta.getAlignedRecordCount();
		        }
		        
		        double[] counts = new double[2*regionWidth+1];
				Arrays.fill(counts, 1.0E6/count*pseudoCount);

				while(sri.hasNext()){
					SAMRecord rec = sri.next();
					int quali=rec.getMappingQuality();
					if(quali>20){
						List<AlignmentBlock> lab = rec.getAlignmentBlocks();
						for(AlignmentBlock ab : lab){
							int blockstart = ab.getReferenceStart();
							int len = ab.getLength();
							for(int i=Math.max(start, blockstart)-start;i<Math.min(end, blockstart+len)-start;i++){
								counts[i]+= 1.0E6/count;
							}
						}
					}
					
				}
				
				if(k<posFiles.length){
					en.addPositives(counts);
				}else{
					en.addNegatives(counts);
				}
				
			}
			
			
		}
		
		LinkedList<Result> ress = new LinkedList<>();
		
		StringBuffer sb = new StringBuffer();
		
		sb.append("#"+Entry.getHeader()+"\tlog fold-change\tcenter-max\n");
		
		for(Entry en : predsList){
			Pair<double[], String> pair = en.getProfileResult(files, windowWidth, t, compare);
			
			double[] ri = pair.getFirstElement();
			if(pair.getSecondElement()!=null){
				sb.append(en+"\t"+ri[0]+"\t"+((int)ri[1])+"\n");
				ress.add( new TextResult("Profile for "+en.chr+":"+en.pos+":"+en.strand, "", new FileParameter.FileRepresentation("", pair.getSecondElement()), "tsv", getToolName(), null, true) );
			}
		}
		
		TextResult tr = new TextResult("Differentially abundant", "", new FileParameter.FileRepresentation("", sb.toString()), "tsv", getToolName(), null, true);
		
		ress.addFirst(tr);
		
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(ress.toArray(new Result[0])), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	private String[] getFiles(ExpandableParameterSet value) {
		String[] vals = new String[value.getNumberOfParameters()];
		for(int i=0;i<value.getNumberOfParameters();i++){
			vals[i] = ((ParameterSet)value.getParameterAt(i).getValue()).getParameterAt(0).getValue().toString();
		}
		return vals;
	}
	
	@Override
	public String getToolName() {
		return "DerTALE";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "dertale";
	}

	@Override
	public String getDescription() {
		return "filters genome-wide predictions for differential expression";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( DerTALE.class.getClassLoader().getResourceAsStream( "projects/tals/rnaseq/DerTALE.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}
	
	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {
		
	}

	@Override
	public String[] getReferences() {
		return new String[] {"@article{erkes19preditale,\n" + 
				"	title = {{PrediTALE}: A novel model learned from quantitative data allows for new perspectives on {TALE} targeting},\n" + 
				"	author = {Erkes, Annett AND M\\\"ucke, Stefanie AND Reschke, Maik AND Boch, Jens AND Grau, Jan},\n" + 
				"	journal = {PLOS Computational Biology},\n" + 
				"	year = {2019},\n" + 
				"	volume = {15},\n" + 
				"	number = {7},\n" + 
				"	pages = {1-31},\n" + 
				"	doi = {10.1371/journal.pcbi.1007206}\n" + 
				"	}\n"};
	}


}
