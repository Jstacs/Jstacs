package projects.tals.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.tals.prediction.QuickTBSPredictionTool;
import umontreal.ssj.probdist.FisherFDist;
import de.jstacs.DataType;
import de.jstacs.algorithms.optimization.ConstantStartDistance;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
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
import de.jstacs.utils.Normalisation;
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

public class DerTALEv2 implements JstacsTool {

	private enum Compare{
		EXTREMES,
		MEDIAN,
		MEAN
	}
	
	private static class Entry{
		private String chr;
		private int pos;
		private String strand;
		private LinkedList<CandidateRegion> candidateRegions;
		private String sequence; // NU: Für das gff3
		private String tal; // NU: Für das gff3
		
		public static String getHeader(){
			return "Chr\tPosition\tStrand";
		}
		
		public String toString(){
			return chr+"\t"+pos+"\t"+strand;
		}
		
		public Entry(String chr, int pos, String strand){
			this.chr = chr;
			this.pos = pos;
			this.strand = strand;
		}
		
		// NU: Für das gff3:
		public Entry(String chr, int pos, String strand, String sequence, String tal){
			this.chr = chr;				
			this.pos = pos;				
			this.strand = strand;		
			this.sequence = sequence;	
			this.tal = tal;				
		}
		
		public String getTal() {		
			return tal;					
		}								
		
		public String getSequence() {	
			return sequence;			
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
		
		public void addCandidateRegion(CandidateRegion candReg){
			if(this.candidateRegions == null){
				this.candidateRegions = new LinkedList<CandidateRegion>();
			}
			this.candidateRegions.add(candReg);
		}
		
		public LinkedList<CandidateRegion> getCandidateRegions() {
			return candidateRegions;
		}

	}
	
	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new DerTALEv2(),new QuickTBSPredictionTool());
		
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
			pars.add(new SimpleParameter(DataType.INT, "Region width", "Number of bases around the predicted site", true,500));
			pars.add(new SimpleParameter(DataType.DOUBLE, "Threshold", "Threshold on the log differential abundance", true,1.0));
			pars.add(new EnumParameter(Stranded.class, "Defines whether the reads are stranded. "
					+ "In case of FR_FIRST_STRAND, the first read of a read pair or the only read in case of single-end data is assumed to be located on forward strand of the cDNA, i.e., reverse to the mRNA orientation. "
					+ "If you are using Illumina TruSeq you should use FR_FIRST_STRAND."
					, true ));
			pars.add(new SimpleParameter(DataType.INT, "Coverage cutoff", "Minimum amount of reads as coverage cuttoff.", true,10));
			pars.add(new SimpleParameter(DataType.INT, "Region elongation value", "Amount of bases a region is elongated if coverage is above half of coverage cuttoff at start/end of region.", true,100));
			pars.add(new SimpleParameter(DataType.INT, "Minimum length of candidate region", "Minimum length of candidate region.", true,100));
			pars.add(new SimpleParameter(DataType.INT, "Minimum coverage of the Candidate Region", "Minimum coverage of the Candidate Region.", true,50));
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
		//System.out.println(predsFile.getFilename());
		
		String[] posFiles = getFiles((ExpandableParameterSet)parameters.getParameterAt(1).getValue());
		String[] negFiles = getFiles((ExpandableParameterSet)parameters.getParameterAt(2).getValue());
		String[] files = new String[posFiles.length+negFiles.length];
		
		
		System.arraycopy(posFiles, 0, files, 0, posFiles.length);
		System.arraycopy(negFiles, 0, files, posFiles.length, negFiles.length);
		
		String[] fileNames = new String[files.length];
		for(int i=0;i<fileNames.length;i++){
			fileNames[i]=Paths.get(files[i]).getFileName().toString();
		}
		
		int numberOfSamples=files.length;		
		
		int numTop = (int) parameters.getParameterAt(3).getValue(); //100
		int regionWidth = (int) parameters.getParameterAt(4).getValue(); //500
		double t = (double) parameters.getParameterAt(5).getValue();
		
		Stranded stranded = (Stranded) parameters.getParameterAt(6).getValue();//TODO anpassen für jedes bam separat angeben?
		int covCuttoff=(int) parameters.getParameterAt(7).getValue(); //10
		int regionElongation=(int) parameters.getParameterAt(8).getValue(); //100
		int minLengthRegion=(int) parameters.getParameterAt(9).getValue(); //100	
		int minCov=(int) parameters.getParameterAt(10).getValue(); //NU: 50 
		
		LinkedList<Entry> predsList = new LinkedList<>(); //einlesen der Liste an Predictions für einen TAL
		BufferedReader br = new BufferedReader(new StringReader(predsFile.getContent())); //BXOR1_RVDs/Predicted_binding_sites_for_TalAD9.tsv
		//		# Seq-ID	Position	Strand	Score	Sequence	Approx. p-value	RVDs	TALE
		//NU:	# Seq-ID	Position	Distance to end	Strand	Score	Sequence	Approx. p-value	RVDs	TALE		(Neuer PrediTALE-Output)
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
				//PrediTALE
			}else if (curr.startsWith("#")) {
				protocol.append("Reading PrediTALE predictions...\n");
				while( (curr = br.readLine()) != null ){
					
					if(!curr.startsWith("#")){
						String[] parts = curr.split("\t");
						predsList.add(new Entry(parts[0], Integer.parseInt(parts[1]), parts[2], parts[4], parts[7]));	// NU: Strand ist jetzt in vierter Spalte (s.o.) und für gff3 wird parts[5], parts[8] zusätzlich eingelesen
						
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
		
		
		SamReader sr[]=new SamReader[numberOfSamples];
		double libSize[]=new double[numberOfSamples];
		Arrays.fill(libSize, 0.0);
		double meanLibSize=0.0;
		
		for( int k = 0; k < numberOfSamples; k++ ) {
			String fName = files[k];
			//System.out.println("k. bam file: "+k);
			
			if(!new File(fName+".bai").exists()){
				throw new RuntimeException("No index found for file "+fName+". The index must be in the same directory as the specified BAM file with filename "+fName+".bai.");
			}
			protocol.append("Opening mapping file "+fName+"...\n");
			//SamReader of k.th samlpe 
			sr[k] = srf.open(SamInputResource.of(new File(fName)).index(new File(fName+".bai")));
		}
		protocol.append("...calculate scaling factor for each library...\n");
		//calculate library size - sum of all mapped reads of all chromosomes 
		for( int k = 0; k < numberOfSamples; k++ ) {
			
			int nRefs=sr[k].getFileHeader().getSequenceDictionary().size();
			for(int i=0;i<nRefs;i++){
				libSize[k]+=sr[k].indexing().getIndex().getMetaData(i).getAlignedRecordCount();
			}
			//System.out.println("libSize "+k+": "+libSize[k]);
			meanLibSize+=libSize[k];
		}
		
		meanLibSize=ToolBox.mean(libSize);
		
		//System.out.println("meanLibSize: "+meanLibSize);
		
	 	//calculate scaling factor for each library
		double[] scalingFactor=new double[numberOfSamples];
		Arrays.fill(scalingFactor, 0.0);
		for( int k = 0; k < numberOfSamples; k++ ) {
			scalingFactor[k]=meanLibSize/libSize[k];
			//System.out.println("scalingFactor "+k+": "+scalingFactor[k]);
		}
		protocol.append("Searching for differentially expressed regions...\n");
		StringBuffer statValues=new StringBuffer();
		statValues.append("Chr"+"\t"+"BS-Pos"+"\t"+"BS strand"+"\t"+"start region"+"\t"+"end region"+"\t"+"strand"+"\t"+"isSignificant"+"\t"+"Fstat"+"\t"+"threshold"+"\t"+"p-value"+"\t"+"lfc"+"\n");
		for(Entry en : predsList){
			//System.out.println("BS-Pos: "+en.getPos());
			TreeMap<Integer,Integer>counts[][]=new TreeMap[2][numberOfSamples]; //[Fwd or Bwd counts][#files] 
			boolean isCovAboveThresholdEnd=false;
			boolean isCovAboveThresholdStart=false;
			
			for( int k = 0; k < numberOfSamples; k++ ) {
				int start = Math.max(0, en.getPos()-regionWidth); //500 positionen davor und 500 Positionen danach
				int end = Math.min(sr[k].getFileHeader().getSequence(en.chr).getSequenceLength(),  en.getPos()+regionWidth );
				
				counts[0][k]=new TreeMap<Integer,Integer>(); 
				counts[1][k]=new TreeMap<Integer,Integer>(); 
				
				calculateCovergage(sr[k], stranded, start, end, en, counts[0][k], counts[1][k]); //calculate coverage at each pos +-regionwidth (500)
				
				
					for(int fwd_bwd=0;fwd_bwd<2;fwd_bwd++){
						if((counts[fwd_bwd][k].get(counts[fwd_bwd][k].firstKey())>(covCuttoff/2))){
							if(counts[0][k].firstKey()>0){ //chr is longer
								isCovAboveThresholdStart=true;
							}
						}
						if((counts[fwd_bwd][k].get(counts[fwd_bwd][k].lastKey())>(covCuttoff/2))){
							if(counts[0][k].lastKey()<sr[k].getFileHeader().getSequence(en.chr).getSequenceLength()){//chr is longer
								isCovAboveThresholdEnd=true;
							}
						}
					}
					
			}
			double halfCuttoff=(double)covCuttoff/2.0;

			while(isCovAboveThresholdEnd){	
				isCovAboveThresholdEnd=false;
				for( int k = 0; k < numberOfSamples; k++ ) {
					int start = Math.max(0, counts[0][k].lastKey()+1);
					int end = Math.min(sr[k].getFileHeader().getSequence(en.chr).getSequenceLength(),  start+regionElongation );
					
					//System.out.println(start);
					//System.out.println(end);
					calculateCovergage(sr[k], stranded, start, end, en, counts[0][k], counts[1][k]);
					
					if((counts[0][k].get(counts[0][k].lastKey())>halfCuttoff)||(counts[1][k].get(counts[1][k].lastKey())>halfCuttoff)){
						if(counts[0][k].lastKey()<sr[k].getFileHeader().getSequence(en.chr).getSequenceLength()){//chr is longer
							isCovAboveThresholdEnd=true;
						}
					}
					
				}
				
			}
		
			while(isCovAboveThresholdStart){
				
				isCovAboveThresholdStart=false;
				for( int k = 0; k < numberOfSamples; k++ ) {
					int end = Math.max(0,  counts[0][k].firstKey()-1 );
					int start = Math.max(0, end-regionElongation);
					
					calculateCovergage(sr[k], stranded, start, end, en, counts[0][k], counts[1][k]);
					
					if((counts[0][k].get(counts[0][k].firstKey())>halfCuttoff)||(counts[1][k].get(counts[1][k].firstKey())>halfCuttoff)){
						if(counts[0][k].firstKey()>0){ //chr is longer
							isCovAboveThresholdStart=true;
						}
						
					}
					
				}
				
			}
			
			ArrayList<CandidateRegion> candRegions[]=new ArrayList[2];
			candRegions[0]=new ArrayList<CandidateRegion>(); //Fwd
			candRegions[1]=new ArrayList<CandidateRegion>();//Bwd
			
			// NU: Es wird nur noch auf dem Strang nach einer Region gesucht, auf dem auch die Target-Box liegt:
			if(en.getStrand().equals("+")) { 
				//System.out.println("...calc candRegionsFwd...");
				getCandidateRegions(candRegions[0],counts[0],posFiles.length,negFiles.length,covCuttoff,minLengthRegion,en.getPos(),regionWidth,Strand.FWD,en.getStrand(),t,minCov);
			}else {
				if(en.getStrand().equals("-")) {
					//System.out.println("...calc candRegionsBwd..."); 
					getCandidateRegions(candRegions[1],counts[1],posFiles.length,negFiles.length,covCuttoff,minLengthRegion,en.getPos(),regionWidth,Strand.REV,en.getStrand(),t,minCov);
				}else {	// NU: wenn Strang unknown --> kann das passieren?
					//System.out.println("...calc candRegionsFwd...");
					getCandidateRegions(candRegions[0],counts[0],posFiles.length,negFiles.length,covCuttoff,minLengthRegion,en.getPos(),regionWidth,Strand.FWD,en.getStrand(),t,minCov);
					
					//System.out.println("...calc candRegionsBwd...");
					getCandidateRegions(candRegions[1],counts[1],posFiles.length,negFiles.length,covCuttoff,minLengthRegion,en.getPos(),regionWidth,Strand.REV,en.getStrand(),t,minCov);
				}
			}
			
			for(int fwd_bwd=0;fwd_bwd<2;fwd_bwd++){
				if(!candRegions[fwd_bwd].isEmpty()){
					for (CandidateRegion candReg : candRegions[fwd_bwd]){
						//System.out.println("start: "+candReg.start+", stop: "+", strand: "+candReg.strand);
						// calcualte reads within Candidate Region FWD
						countReadsWithinRegion(sr, posFiles.length,negFiles.length, candReg, en,stranded);
						double[] yCounts=new double[numberOfSamples];
						Arrays.fill(yCounts, 0.0);
						int[][] xIndicators=new int[numberOfSamples][2]; // two condictions (treatment vs mock)
						for( int k = 0; k < numberOfSamples; k++ ) {
							Arrays.fill(xIndicators[k],0);
						}
						
						for( int k = 0; k < numberOfSamples; k++ ) {
						
							candReg.setNormalizedCounts(k, (candReg.getCountReads(k)+1.0)*scalingFactor[k]);				
							yCounts[k]=Math.log(candReg.getNormalizedCount(k));												

							if(k<posFiles.length){
								xIndicators[k][1]=1;
								//System.out.println("Treatment");
								//System.out.println(candReg.getCountReads(k)+"->"+candReg.getNormalizedCount(k));
								
							}else{
								xIndicators[k][0]=1;
								//System.out.println("Control");
								//System.out.println(candReg.getCountReads(k)+"->"+candReg.getNormalizedCount(k));
							}
						}
						double lfc = Normalisation.getLogSum(0, posFiles.length, yCounts) - Normalisation.getLogSum(posFiles.length, yCounts.length, yCounts);
						//calculate RSS1
						DerTALEv2ComplexFunction df=new DerTALEv2ComplexFunction(yCounts, xIndicators);
						double meanCounts=ToolBox.mean(yCounts);
						//System.out.println("meanCounts: "+meanCounts);
						double[] currentParameters={meanCounts,0.0,0.0};
						double linEps=1E-6;
						Optimizer.optimize( Optimizer.CONJUGATE_GRADIENTS_PRP, df, currentParameters, new SmallDifferenceOfFunctionEvaluationsCondition( 1E-6 ), linEps, new ConstantStartDistance( 1E-4 ), System.out );
						//System.out.println(Arrays.toString(currentParameters));
						//calculate F
						double RSS1=df.evaluateFunction(currentParameters);
						//System.out.println(RSS1);
						
						//calculate RSS0
						double RSS0=0.0;
						for( int k = 0; k < numberOfSamples; k++ ) {
							double temp=yCounts[k]-meanCounts;
							RSS0+=temp*temp;
						}
						//System.out.println("RSS0: "+RSS0);
						double offset=0.0;
						//System.out.println("                                                                                               xIndicators[0].length:"+xIndicators[0].length);
						double Fstat=((RSS0-RSS1)/(xIndicators[0].length-1))/(offset+(RSS1/(numberOfSamples-xIndicators[0].length)));
						//System.out.println("Fstat:"+Fstat);
						FisherFDist f = new FisherFDist(xIndicators[0].length-1, numberOfSamples-xIndicators[0].length);
						double pval = f.barF(Fstat);
						double threshold = f.inverseF(1-0.05); //1-alpha //TODO signifikanzniveau alpha als parameter!
						//System.out.println("pval: "+pval);
						candReg.setpValue(pval);
						//System.out.println("threshold: "+threshold);
						String strand="-";
						if(fwd_bwd==0){
							strand="+";
						}
						boolean isSignificant=false;
						if(Fstat>threshold){
							isSignificant=true;
							//System.out.println("significant: BS:"+en.getChr()+":"+en.getPos()+":"+en.strand+", candRegion: "+candReg.getStart()+":"+candReg.getEnd()+":"+strand );
						}/*else{
							System.out.println("not significant! BS:"+en.getChr()+":"+en.getPos()+":"+en.strand+", candRegion: "+candReg.getStart()+":"+candReg.getEnd()+":"+strand );
						}*/
						candReg.setSignificance(isSignificant);
						statValues.append(en.getChr()+"\t"+en.getPos()+"\t"+en.strand+"\t"+candReg.getStart()+"\t"+candReg.getEnd()+"\t"+strand+"\t"+isSignificant+"\t"+Fstat+"\t"+threshold+"\t"+pval+"\t"+lfc+"\n");
						
						if(Fstat>threshold){
							for( int k = 0; k < files.length; k++ ) {
								String fName = files[k];
								
								//System.out.println("...calculateProfileCounts...");
								double[] profileCounts=calculateProfileCounts(en, srf, fName, candReg, scalingFactor[k], stranded);
								
								//System.out.println(Arrays.toString(profileCounts));
								if(k<posFiles.length){
									candReg.addPositives(profileCounts);
								}else{
									candReg.addNegatives(profileCounts);
								}
								
							}
							
							en.addCandidateRegion(candReg);
							
						}
							
						
							
						}
					}
				}
			
		}
		protocol.append("Finished.\nPreparing output...");
		
		LinkedList<Result> ress = new LinkedList<>();
		
		ArrayList <String> results = new ArrayList<>();	//NU
		String gff3Content = "##gff-version 3";	// NU
		
		for(Entry en : predsList){
			int countCands=0;
			if(en.getCandidateRegions()!=null){
				//System.out.println("...getProfileResult for significant regions...");
				for(CandidateRegion candReg : en.getCandidateRegions()){
					countCands++;
					String profileResult = candReg.getProfileResult(fileNames, t, en.getPos(),candReg);
					
					//NU: Wenn schon ein Profil für die gleiche Target-Box mit der gleichen Region erstellt wurde, dann wird nicht noch ein zweites dafür erstellt (verhindert Fehlermeldung für doppelte Profile durch aberrante Repeats)
					if(profileResult!=null && results.contains(en.chr+":"+en.pos+":"+en.strand+":"+candReg.getStart()+":"+candReg.getEnd()+":"+candReg.getStrand().toString()) == false){ 
						ress.add( new TextResult("Profile for "+en.chr+":"+en.pos+":"+en.strand+" candidateRegion:"+candReg.getStart()+":"+candReg.getEnd()+":"+candReg.getStrand().toString(), "", new FileParameter.FileRepresentation("", profileResult), "tsv", getToolName(), null, true) );	
						results.add(en.chr+":"+en.pos+":"+en.strand+":"+candReg.getStart()+":"+candReg.getEnd()+":"+candReg.getStrand().toString()); 	//NU
						
						// NU: Daten der Region wird an gff3-String angehangen:
						gff3Content += "\n"+en.chr+"\t"+getToolName()+"\t.\t"+candReg.getStart()+"\t"+candReg.getEnd()+"\t.\t"+candReg.getStrand().toString()+"\t.\tName="+en.pos+"_candidateRegion:"+candReg.getStart()+"_"+candReg.getEnd();
						
						// NU: Daten der Target-Box ist nur für PrediTALE-Eingabe im gff3:
						if(en.getSequence() != null && en.getTal() != null) {	
							gff3Content += "\n"+en.chr+"\t"+getToolName()+"\t.\t"+(en.pos+1)+"\t"+(en.pos+en.getSequence().length())+"\t.\t"+en.strand+"\t.\tName="+en.getTal()+"_"+en.pos;
						}
					}
				}
			}
			//System.out.println("Found "+countCands+" significant regions for BS-Pos: "+en.pos);
		}
		
		
		// NU: Erstellen und Hinzufügen der gff3 zur Ausgabe:
		TextResult gff3T =  new TextResult("Tals_and_regions", "", new FileParameter.FileRepresentation("", gff3Content), "gff3", getToolName(), null, true);
		ress.addFirst(gff3T);
		
		TextResult trstatValues = new TextResult("Differentially abundant", "", new FileParameter.FileRepresentation("", statValues.toString()), "tsv", getToolName(), null, true);
		
		ress.addFirst(trstatValues);
		protocol.append("Finished.\nWriting output...\n");
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(ress.toArray(new Result[0])), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}
	
	private boolean isFwd(SAMRecord sr, Stranded stranded) {
		if(stranded == Stranded.FR_UNSTRANDED) {
			return true;
		}else if(stranded == Stranded.FR_FIRST_STRAND){
			if( 
					(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag() ) ||
					(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag() ) ||
					( (!sr.getReadPairedFlag() || !sr.getProperPairFlag()) && !sr.getReadNegativeStrandFlag() )
				) {
				return true;
			}else {
				return false;
			}
		}else {
			if( 
					(sr.getReadPairedFlag() && sr.getFirstOfPairFlag() && sr.getReadNegativeStrandFlag() ) ||
					(sr.getReadPairedFlag() && sr.getSecondOfPairFlag() && sr.getMateNegativeStrandFlag() ) ||
					( (!sr.getReadPairedFlag() || !sr.getProperPairFlag()) && sr.getReadNegativeStrandFlag() )
				) {
				return false;
			}else {
				return true;
			}
		}
	}
	
	private void calculateCovergage(SamReader sr, Stranded stranded, int start, int end, Entry en, TreeMap<Integer,Integer> countsFwd, TreeMap<Integer,Integer> countsBwd) {
		for(int i=start;i<=end;i++){
			countsFwd.put(i, 0);
			countsBwd.put(i, 0);
			
		}
		
		SAMRecordIterator sri;
		sri = sr.query(en.chr, start, end, false);	// NU: false: nur overlapping mit mind. eine Pos., true: vollständig enthalten; false verursacht teilweise Probleme bei der Laufzeit
		boolean isReadFwd;
		int readGenomicStart=-1;
		int readGenomicEnd=-1;
		
		while(sri.hasNext()){
		
			SAMRecord rec = sri.next();
				
			int quali=rec.getMappingQuality();
			
			if(quali>20){
				System.out.flush();
				isReadFwd=isFwd(rec, stranded); // is read forward
				
				readGenomicStart=rec.getAlignmentStart();
				readGenomicEnd=rec.getAlignmentEnd();
				
				int readStartRegion = Math.max(readGenomicStart, start); //500 positionen davor und 500 Positionen danach
				int readEndRegion = Math.min(readGenomicEnd,  end );
				
				if(isReadFwd){ 
					for(int i=readStartRegion;i<=readEndRegion;i++){
						countsFwd.replace(i, countsFwd.get(i)+1);
					}
				}else{
					for(int i=readStartRegion;i<=readEndRegion;i++){
						countsBwd.replace(i, countsBwd.get(i)+1);
					}
				}
			}
		}
		sri.close();
	}
	
	private void countReadsWithinRegion(SamReader[] sr, int posFilesLength, int negFilesLength,CandidateRegion candReg, Entry en,Stranded stranded){
		boolean isCandRegionFWD=true;
		if(candReg.getStrand()==Strand.REV){
			isCandRegionFWD=false;
		}
		
		for( int k = 0; k < sr.length; k++ ) {
			SAMRecordIterator sri;
			sri = sr[k].query(en.chr, candReg.start, candReg.end, false);
			
			boolean isReadFwd;

			while(sri.hasNext()){
				SAMRecord rec = sri.next();
				
				int quali=rec.getMappingQuality();
				if(quali>20){
					
					isReadFwd=isFwd(rec, stranded); // is read forward
					if(isReadFwd==isCandRegionFWD){
						candReg.setCountReads(k, candReg.getCountReads(k)+1);
					}
				}
			}
			sri.close();
		}
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
	
	private void getCandidateRegions(ArrayList<CandidateRegion> candRegions,TreeMap<Integer,Integer>[]counts,int posFilesLength, int negFilesLength, int covCuttoff, int minLengthRegion,int BSpos,int regionWidth, Strand strand, String TStrand, double t, int minCov){     
		// NU: strand --> Strang der Region, TStrand --> Strang der Target-Box
		int candStart=-1;
		int candEnd=-1;
		boolean isWithinCand=false;
		/*System.out.println("BSpos: "+BSpos);
		System.out.println("FistPos: "+counts[0].firstKey());
		System.out.println("LastPos: "+counts[0].lastKey());
		System.out.println("keySet-size: "+counts[0].keySet().size());
		System.out.println("keySet: "+counts[0].keySet().toString());*/
		
		double mean30PositionsPos = 0;	// NU
		double mean30PositionsNeg = 0;	// NU
		double max = 0;					// NU
		double BSposCount = 0;			// NU
		int sumAboveThreshold = 0;		// NU
		
		for(int pos : counts[0].keySet()){
			double meanPos=0;
			double meanNeg=0;
			for( int k = 0; k < counts.length; k++ ) {
				
				if(k<posFilesLength){ //treatment
					meanPos+=(double)counts[k].get(pos);							
				
				}else{//control
					meanNeg+=(double)counts[k].get(pos);							
				
				}
			}
			
			meanPos=meanPos/posFilesLength;											
			meanNeg=meanNeg/negFilesLength;
			
			if(((meanPos>=covCuttoff)||(meanNeg>=covCuttoff))){																					

				if(isWithinCand==false){											
					candStart=pos;													
					if(candStart>BSpos){//TAL-BS is before candStart							
						if(Math.abs(pos-BSpos)>=regionWidth){						
							isWithinCand=false;										
							break;															
						}
					}
					
				}
				candEnd=pos;
				isWithinCand=true;
				
				// NU: Maximum der Region speichern:
				if(meanPos > max) {												
					max = meanPos;
				}															
				if(meanNeg > max) {													
					max = meanNeg;
				}
				
				// NU: Speichern, wie oft LfC über dem threshold:
				if(Math.abs(Math.log(meanPos+1) - Math.log(meanNeg+1)) > t) {		
					sumAboveThreshold++;
				}
				
			}else{ // 															
				if(isWithinCand){													
					if(candEnd<BSpos){ //TAL-BS after candEnd						
						if(Math.abs(BSpos-pos)>=regionWidth){						
							isWithinCand=false;										
							mean30PositionsPos = 0; 								// NU				
							mean30PositionsNeg = 0;									// NU							
							max = 0; 												// NU							
							sumAboveThreshold = 0;									// NU								
							continue; 												// NU vorher: break;												
						}															
					}										
			
					if((candEnd-candStart)>=minLengthRegion){  						
						
						// NU: Berechnung des Mittelwerts der 30Pos vor der Target-Box, wenn die Target-Box auf dem - Strang:
						if(TStrand.equals("-")) {									
							for(int i = BSpos-29; i <= BSpos; i++) {
								meanPos=0;
								meanNeg=0;
								for( int k = 0; k < counts.length; k++ ) {
									if(k<posFilesLength){ //treatment
										meanPos+=(double)counts[k].get(i);			// NU: aufaddieren aller counts für treatment aus jeweils den drei bam an der Poition pos 
									}else{//control
										meanNeg+=(double)counts[k].get(i);			// NU: aufaddieren aller counts für control aus jeweils den drei bam an der Poition pos 			
									}
								}
								
								mean30PositionsPos += meanPos/posFilesLength;// treatment
								mean30PositionsNeg += meanNeg/negFilesLength;// control
								
								// NU: Speichern der (größeren) Coverage an der Position der Target-Box:
								if(i == BSpos) {									
									if(meanPos/posFilesLength > meanNeg/negFilesLength) {
										BSposCount = meanPos/posFilesLength;			
									}else {
										BSposCount = meanNeg/negFilesLength;
									}
								}
							}
						}else {														
							// NU: Berechnung des Mittelwerts der 30Pos nach der Target-Box, wenn die Target-Box auf dem + Strang:
							for(int i = BSpos; i < BSpos+30; i++) {
								meanPos=0;
								meanNeg=0;
								for( int k = 0; k < counts.length; k++ ) {
									if(k<posFilesLength){ //treatment
										meanPos+=(double)counts[k].get(i);			// NU: aufaddieren aller counts für treatment aus jeweils den drei bam an der Poition pos 
									}else{//control
										meanNeg+=(double)counts[k].get(i);			// NU: aufaddieren aller counts für control aus jeweils den drei bam an der Poition pos 			
									}
								}
								
								mean30PositionsPos += meanPos/posFilesLength;// treatment
								mean30PositionsNeg += meanNeg/negFilesLength;// control
								
								// NU: Speichern der (größeren) Coverage an der Position der Target-Box:
								if(i == BSpos) {									
									if(meanPos/posFilesLength > meanNeg/negFilesLength) {
										BSposCount = meanPos/posFilesLength;			
									}else {
										BSposCount = meanNeg/negFilesLength;
									}
								}
							}
						}
						
						// NU: Region muss komplett vor/nach Target-Box liegen   oder   30Pos-Fenster vor/nach Target-Box liegt unter cuttoff
						if((TStrand.equals("+") && (BSpos+30) < candStart) || (TStrand.equals("-") && (BSpos-30) > candEnd) 
							|| (TStrand.equals("+") && (BSpos+30) >= candStart && candEnd > (BSpos+30) && (mean30PositionsPos/30.0 < covCuttoff && mean30PositionsNeg/30.0 < covCuttoff))
							|| (TStrand.equals("-") && (BSpos-30) <= candEnd && candStart < (BSpos-30) && (mean30PositionsPos/30.0 < covCuttoff && mean30PositionsNeg/30.0 < covCuttoff))) {		
							
							if(BSposCount < (max/2.0) && max > minCov) { 			// NU: An Target-Box muss Coverage kleiner als die Hälfte des maximalen der Region sein und Maximum des Fensters muss größer 50 (default) sein
								if(sumAboveThreshold > ((candEnd-candStart)/5)) {	// NU: LogFoldChange muss an ein Fünftel der Pos in der Region über dem threshold liegen
									//System.out.println("candStart:"+candStart);
									//System.out.println("candEnd:"+candEnd);
									//TODO ende nicht zu weit von TAL entfernt, falls TAL hinter cand liegt!!!
									candRegions.add(new CandidateRegion(candStart, candEnd, strand,posFilesLength,negFilesLength));
								}
							}
						}
					}
				}
				isWithinCand=false;
				mean30PositionsPos = 0;	// NU
				mean30PositionsNeg = 0;	// NU
				max = 0;				// NU
				sumAboveThreshold = 0;	// NU
			}	
			
		}
	}
	
	private static class CandidateRegion{
		private int start;
		private int end;
		private Strand strand;
		private int[] countReads;
		private double[] normalizedCounts;
		private LinkedList<double[]> positives;
		private LinkedList<double[]> negatives;
		private int correctedCandStartPos; //0-FWD candidate, 1-BWD candidate
		private double pValue;
		private boolean isSignificant;

		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}

		public Strand getStrand() {
			return strand;
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

		public int getcorrectedCandStartPos() {
			return correctedCandStartPos;
		}

		public void setcorrectedCandStartPos(int correctedCandStartPos) {
			
			this.correctedCandStartPos = correctedCandStartPos;
		}
		
		public double getpValue(){
			return this.pValue;
		}
		
		public void setpValue(double pValue){
			this.pValue=pValue;
		}
		public void setSignificance(boolean isSignificant){
			this.isSignificant=isSignificant;
		}
		
		private CandidateRegion(int start, int end, Strand strand, int numberTreatment, int numberControl){
			this.start = start;
			this.end = end;
			this.strand = strand;
			this.countReads=new int[numberTreatment+numberControl];
			Arrays.fill(this.countReads, 0);
			this.normalizedCounts=new double[numberTreatment+numberControl];
			Arrays.fill(this.normalizedCounts, 0.0);
			this.correctedCandStartPos=-1;
			this.pValue=-1.0;
			this.isSignificant=false;
			
		}

		public int getCountReads(int index) {
			return countReads[index];
		}

		public void setCountReads(int index, int count) {
			this.countReads[index] = count;
		}

		public double getNormalizedCount(int index) { //index of the sample
			return normalizedCounts[index];
		}

		public void setNormalizedCounts(int index, double normalizedCount) {
			this.normalizedCounts[index] = normalizedCount;
		}
		
		public String getProfileResult(String[] header, double t, int enPos, CandidateRegion candReg) {
			
			int mid = candReg.positives.getFirst().length/2; //3000
			
			StringBuffer sb = new StringBuffer();
			sb.append("Position");
			for(int i=0;i<header.length;i++){
				sb.append("\t"+header[i]);
				
			}
			sb.append("\tabove threshold\t");
			sb.append("pValue:"+candReg.getpValue()+"\n");
			//System.out.println(correctedCandStartPos);
			int aktPos=0;
			for(int i=0;i<candReg.positives.get(0).length;i++){
				aktPos=i+candReg.getcorrectedCandStartPos();
			
				sb.append(aktPos);
	
				for(int j=0;j<candReg.positives.size();j++){
					sb.append("\t"+candReg.positives.get(j)[i]);
				}
				for(int j=0;j<candReg.negatives.size();j++){
					sb.append("\t"+candReg.negatives.get(j)[i]);
				}
				if((aktPos >= candReg.getStart())&&(aktPos <= candReg.getEnd())){   
					sb.append("\tTRUE\n");
				}else{
					sb.append("\tFALSE\n");
				}
			}
			
			/*System.out.println("BS-pos: "+enPos);
			System.out.println("candReg.getStart(): "+candReg.getStart());
			System.out.println("candReg.getEnd(): "+candReg.getEnd());
			System.out.println("lastPos:"+aktPos);
			System.out.println(candReg.getcorrectedCandStartPos());*/

				
			return sb.toString();
				
		}
	}
	
	private double[] calculateProfileCounts(Entry en,SamReaderFactory srf,String fName, CandidateRegion candReg,double scalingFactor, Stranded stranded){
		//todo for all candRegs des entrys for schleife drum rum!
		
		if(!new File(fName+".bai").exists()){
			throw new RuntimeException("No index found for file "+fName+". The index must be in the same directory as the specified BAM file with filename "+fName+".bai.");
		}
		boolean isCandRegionFWD=true;
		if(candReg.getStrand()==Strand.REV){
			isCandRegionFWD=false;
		}
		int start=candReg.getStart();
		int end=candReg.getEnd();
		int extraBases=100;
		//System.out.println("candReg.getStart()"+candReg.getStart());
		if(start>en.getPos()){ //if BS-pos is before candidate
			start=start-(start-en.getPos());
		}
		else if(end<en.getPos()){//if BS-pos is after candidate
			end=end+(en.getPos()-end);
		}
		start=start-extraBases;
		end=end+extraBases;
		
		candReg.setcorrectedCandStartPos(start);
		//System.out.println("BS-pos: "+en.getPos());
		//System.out.println("correctedStart: "+start);
		
		
		SamReader sr = srf.open(SamInputResource.of(new File(fName)).index(new File(fName+".bai")));
		
		SAMRecordIterator sri = sr.query(en.chr, start, end, false);
       
        double[] counts = new double[end-start+1];
		Arrays.fill(counts, 0.0);
		boolean isReadFwd;
		while(sri.hasNext()){
			SAMRecord rec = sri.next();
			int quali=rec.getMappingQuality();
			
			if(quali>20){
				isReadFwd=isFwd(rec, stranded); // is read forward
				if(isReadFwd==isCandRegionFWD){
					List<AlignmentBlock> lab = rec.getAlignmentBlocks();
					for(AlignmentBlock ab : lab){
						int blockstart = ab.getReferenceStart();
						int len = ab.getLength();
						for(int i=Math.max(start, blockstart)-start;i<Math.min(end, blockstart+len)-start;i++){
							counts[i]+= scalingFactor;
						}
					}
				}
				
			}
			
		}
		
		return counts;
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
		return "DerTALEv2";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "dertalev2";
	}

	@Override
	public String getDescription() {
		return "filters genome-wide predictions for differential expression";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( DerTALEv2.class.getClassLoader().getResourceAsStream( "projects/tals/rnaseq/DerTALEv2.txt" ) ).toString();
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
