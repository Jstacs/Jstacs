package projects.tals.epigenetic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DiscreteSequenceEnumerator;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.QuickScanningSequenceScore;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.LargeSequenceReader;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import htsjdk.samtools.util.RuntimeEOFException;
import projects.tals.RVDSequence;
import projects.tals.prediction.QuickTBSPredictionTool;
import umontreal.iro.lecuyer.probdist.NormalDist;

public class QuickTBSPredictionToolMethylAccessibilityAnnotation_fai implements JstacsTool {

	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new QuickTBSPredictionToolMethylAccessibilityAnnotation_fai());
		cli.run(args);
	}

	public QuickTBSPredictionToolMethylAccessibilityAnnotation_fai() {

	}

	@Override
	public ToolParameterSet getToolParameters() {
		//FileParameter model = new FileParameter("Model", "The model (in XML format)", "xml", true);

		FileParameter genome = new FileParameter("Sequences","The sequences (e.g., a genome) to scan for binding sites","fa,fas,fasta",true);

		FileParameter background = new FileParameter("Background sequences","The sequences (e.g., a genome) for determining the prediction threshold","fa,fas,fasta",true);

		try {

			SelectionParameter bg = new SelectionParameter(DataType.PARAMETERSET,new String[]{"sub-sample","background sequences"}, new ParameterSet[]{
					new SimpleParameterSet(),
					new SimpleParameterSet(background)
			},"Background sample","The sequences for determining the prediction threshold. Either a sub-sample of the input sequences or a dedicated background data set.",true);

			SimpleParameter pvalue = new SimpleParameter(DataType.DOUBLE, "Significance level", "The significance level for determining the prediction threshold", true, new NumberValidator<Comparable<Double>>(0.0, 1E-2), 1E-4);

			SimpleParameter num = new SimpleParameter(DataType.INT,"Number of sites", "The number of expected binding sites for determining the prediction threshold", true, new NumberValidator<Comparable<Integer>>(1, 1000000),10000);

			SelectionParameter tsel = new SelectionParameter(DataType.PARAMETERSET, new String[]{"significance level","number of sites"}, new ParameterSet[]{
					new SimpleParameterSet(pvalue),
					new SimpleParameterSet(num)
			}, "Threshold specification", "The way of defining the prediction threshold. Either by explicitly defining a significance level or by specifying the number of expected sites", true);


			FileParameter tals = new FileParameter("TALEs", "The RVD sequences of the TALE, separated by dashes, in FastA format","fasta,fas,fa", true);

			
			SelectionParameter strand = new SelectionParameter(DataType.PARAMETERSET, new String[]{"both strands","forward strand","reverse strand"}, new ParameterSet[]{
					new SimpleParameterSet( new SimpleParameter(DataType.DOUBLE, "Reverse penalty", "Penalty for predictions on the reverse strand", true, new NumberValidator<Double>(0.0, Double.MAX_VALUE), 0.01) ),
					new SimpleParameterSet(),
					new SimpleParameterSet()
			}, "Strand", "Prediction target sites on both strands, or the forward or reverse strand", true);
			
			//bismark output.cov.gz
			//<chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated> <count unmethylated>
			//chromosome04    887     887     0       0       4
			//chromosome04    892     892     14.2857142857143        1       6
			//chromosome04    893     893     28.5714285714286        2       5

			FileParameter bismark = new FileParameter("Bismark bedGraph output","The bedGraph output of bismark (file.cov.gz) containig <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>","cov.gz",false);
			
			FileParameter fai = new FileParameter("Genome fasta index file","The fasta index file (.fai) of the genome","fai",true);
			
			//narrowPeakFile
			FileParameter narrowPeak = new FileParameter("NarrowPeak File","The output of a peak caller (all.peaks.narrowPeak)","narrowPeak",false);
			
			//normalized CoverageFile/Pileup-Output
			FileParameter coveragePileup =new FileParameter("normalized pileup output","The normalized output of pileup with values larger than zero (file.txt) containig <chromosome> <position> <coverage>","txt",false);
			
			SimpleParameter cov_before = new SimpleParameter(DataType.INT,"Coverage before value", "Number of positions before binding site in coverage profile", false, new NumberValidator<Comparable<Integer>>(1, 500),300);
			SimpleParameter cov_after = new SimpleParameter(DataType.INT,"Coverage after value", "Number of positions after binding site in coverage profile", false, new NumberValidator<Comparable<Integer>>(1, 500),50);
			
			SelectionParameter calculateCoverageArea = new SelectionParameter(DataType.PARAMETERSET, new String[]{"surround target site","on complete sequence"}, new ParameterSet[]{
					new SimpleParameterSet( cov_before,cov_after),
					new SimpleParameterSet()
			}, "calculate coverage area", "Calculate coverage area surround target site, or on complete sequence", false);
			

			SimpleParameter peak_before = new SimpleParameter(DataType.INT,"Peak before value", "Number of positions before binding site in peak profile", false, new NumberValidator<Comparable<Integer>>(1, 500),300);
			SimpleParameter peak_after = new SimpleParameter(DataType.INT,"Peak after value", "Number of positions after binding site in peak profile", false, new NumberValidator<Comparable<Integer>>(1, 500),50);
			
			
			return new ToolParameterSet(getShortName(),genome,bg,tsel,tals,strand,bismark,fai,narrowPeak,coveragePileup,calculateCoverageArea,peak_before,peak_after);

		} catch (ParameterException e) {
			throw new RuntimeEOFException(e);
		}
	}
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {

		progress.setLast(1.0);
		progress.setCurrent(0.0);

		boolean backgroundSet = ((SelectionParameter)parameters.getParameterAt(1)).getSelected()==1;
		String backgroundPath = backgroundSet ? ( (ParameterSet)parameters.getParameterAt(1).getValue()).getParameterAt(0).getValue().toString() : parameters.getParameterAt(0).getValue().toString();
		String genomePath = parameters.getParameterAt(0).getValue().toString();
		String faiPath = parameters.getParameterAt(6).getValue().toString();
		String bismarkPath = null;
		if(parameters.getParameterAt(5).getValue()!=null){
			bismarkPath = parameters.getParameterAt(5).getValue().toString();
		}
		boolean useMethylationData=false;
		MethylationprofilHashMap methylationProfiles=null;
		if(bismarkPath==null){
			bismarkPath="projects/tals/epigenetic/empty.cov.gz";
			
			methylationProfiles=new MethylationprofilHashMap(faiPath,bismarkPath, 0.0f,0.0f);//withoutMethyl
		}else{
			useMethylationData=true;
			methylationProfiles=new MethylationprofilHashMap(faiPath, bismarkPath, 1.0f,0.0f);//Methyl
		}
		String peakPath=null;
		if(parameters.getParameterAt(7).getValue()!=null){
			peakPath = parameters.getParameterAt(7).getValue().toString();//narrowPeak-File
		}
		
		String coveragePath=null;
		if(parameters.getParameterAt(8).getValue()!=null){
			coveragePath = parameters.getParameterAt(8).getValue().toString();
		}
		int peak_before=(Integer) parameters.getParameterAt(10).getValue();
		int peak_after=(Integer) parameters.getParameterAt(11).getValue();
		
		int startStrand = 0, endStrand = 2;
		double strandPenaltyPerc = 0.0;
	
		if(((SelectionParameter)parameters.getParameterAt(4)).getSelected() == 1){
			endStrand = 1;
		}else if(((SelectionParameter)parameters.getParameterAt(4)).getSelected() == 2){
			startStrand = 1;
		}else{
			strandPenaltyPerc = (double) ((ParameterSet)parameters.getParameterAt(4).getValue()).getParameterAt(0).getValue();
		}
		
		boolean calculateAlwaysOnCompleteSeq=false;
		int cov_before=0;
		int cov_after=0;
		boolean bsb = ((SelectionParameter)parameters.getParameterAt(9)).getSelected()==0; 
		if(bsb){//surround box
			cov_before=(int)((ParameterSet)parameters.getParameterAt(9).getValue()).getParameterAt(0).getValue();
			cov_after=(int)((ParameterSet)parameters.getParameterAt(9).getValue()).getParameterAt(1).getValue();
		}else{//on complete promotor
			calculateAlwaysOnCompleteSeq=true;
		}
		String[][] tals = readTALs(((FileParameter)parameters.getParameterAt(3)).getFileContents());

		LinkedList<Result> talRess = new LinkedList<>();
		
		double fac = 1.0/tals.length;
		double last = 0.0;
		NarrowpeakprofilHashMap narrowPeakProfiles=null;
		if(peakPath!=null){
			narrowPeakProfiles=new NarrowpeakprofilHashMap(faiPath,peakPath, peak_before,peak_after,"genome");
		}
		
		//System.out.println("START!");
		//Instant start = Instant.now();
		PileupCoverageprofilHashMap pileupCoverageProfiles=null;
		if(coveragePath!=null){
			pileupCoverageProfiles=new PileupCoverageprofilHashMap(faiPath,coveragePath, cov_before,cov_after,calculateAlwaysOnCompleteSeq);
		}
		
		//Instant end = Instant.now();
		//System.out.println("END: "+Duration.between(start, end));
		
		for(int ta=0;ta<tals.length;ta++){

			String rvdStr = tals[ta][0];
			String talName = tals[ta][1];

			String[] rvds = rvdStr.split("-");
			IntList idxs = new IntList();
			for(int i=0;i<rvds.length;i++){
				if(!rvds[i].toUpperCase().equals(rvds[i])){
					idxs.add(i);
				}
			}

			protocol.appendHeading("Starting predictions for "+talName+"...\n");
			protocol.append("Using "+(backgroundSet? " background set.\n" : " sub-sample of input data.\n"));

			boolean byp = ((SelectionParameter)parameters.getParameterAt(2)).getSelected()==0;
			double p_value = 0;
			if(byp){
				p_value = (Double) ((ParameterSet)parameters.getParameterAt(2).getValue()).getParameterAt(0).getValue();
			}else{
				int nsites = (Integer) ((ParameterSet)parameters.getParameterAt(2).getValue()).getParameterAt(0).getValue();
				p_value = nsites/(double)(new File(genomePath)).length()/2.0;
			}

			if(idxs.length()>0){
				protocol.append("Found aberrant repeats. Correcting p-value for multiple testing ("+Math.pow(2.0, idxs.length())+")\n\n");
			}
			p_value /= Math.pow(2.0, idxs.length());

			protocol.append("Significance level: "+p_value+"\n");
			
			double subsamp = 0;
			if(backgroundSet){
				subsamp = 1;
			}else{
				subsamp = 1E6/(double)(new File(backgroundPath)).length();
			}

			AlphabetContainer alphabet12 = new AlphabetContainer(new DiscreteAlphabet(false, "A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"));
			AlphabetContainer alphabet13 = new AlphabetContainer(new DiscreteAlphabet(false, "A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*"));

			String[] sepRVDs = new String[]{"HD", "NN", "NG", "NI"};
			
			DiscreteSequenceEnumerator en = new DiscreteSequenceEnumerator(new AlphabetContainer(new DiscreteAlphabet(true, "y", "n")), idxs.length(), false);

			double fac2 = fac / Math.pow(2.0, idxs.length());

			LinkedList<ComparableElement<ResultSet,Double>> ll = new LinkedList<>();
			
			while(en.hasMoreElements()){
				Sequence temp = en.nextElement();
				String[] curr = rvds.clone();
				for(int i=0;i<idxs.length();i++){
					if(temp.discreteVal(i)==0){
						curr[idxs.get(i)] = "";
					}
				}
				String currS = String.join("-", curr).replaceAll("-+", "-");

				RVDSequence eff = new RVDSequence(alphabet12, alphabet13, currS);

				LFModularConditional9CExtMethyl lfmod = new LFModularConditional9CExtMethyl( FileManager.readFile(  "projects/tals/prediction/preditale_quantitative_PBM.xml"), alphabet13, RVDSequence.getContainerRVD(alphabet12, alphabet13), sepRVDs);
				
				double[][] pwm = lfmod.toPWM(eff);
				
				double max = 0;
				double min = 0;
				for(int i=0;i<pwm.length;i++){
					max += ToolBox.max(0,4,pwm[i]);
					min += ToolBox.min(0,4,pwm[i]);
				}
				//protocol.append("min: "+min+", max: "+max+", diff: "+(max-min)+"\n");
				double diff = (max-min);
				double strandPenalty = -diff*strandPenaltyPerc;
				protocol.append("Effective strand penalty: "+strandPenalty+"\n");

				QuickScanningSequenceScore model = new PFMWrapperTrainSMMethyl(DNAAlphabetContainer.SINGLETON, "", pwm);

				int kmer = Math.min(10, model.getLength()*2/3);
				int[] starts = new int[]{0,(model.getLength()-kmer)/3,(model.getLength()-kmer)*2/3,model.getLength()-kmer};
				protocol.append("Target site length: "+model.getLength()+"\n");
				protocol.append("Using "+kmer+"-mers starting at positions "+Arrays.toString(starts)+"\n");

				NormalDist nd = getThreshold(backgroundPath,model,subsamp,startStrand,endStrand,strandPenalty); 
				progress.setCurrent(last + 0.3*fac2);

				double t = nd.inverseF(1.0-p_value);

				protocol.append("Effective threshold: "+t+"\n");

				boolean[][] use = model.getInfixFilter(kmer, t, starts);

				double[] us = new double[use.length];

				for(int i=0;i<use.length;i++){
					for(int j=0;j<use[i].length;j++){
						if(use[i][j]){
							us[i]++;
						}
					}
				}
				
				protocol.append("Number of "+kmer+"-mers passing filtering: "+Arrays.toString(us)+"\n\n");
				
				int[] o = ToolBox.order(us, false);
				boolean[][] temp2 = new boolean[use.length][];
				int[] temps = new int[starts.length];
				for(int i=0;i<use.length;i++){
					temp2[i] = use[o[i]];
					temps[i] = starts[o[i]];
				}
				use = temp2;
				starts = temps;

				protocol.appendHeading("Predicting sites for RVD sequence "+eff.toString("-", 0, eff.getLength())+"...\n");
				getSites(ll, eff, talName, useMethylationData,methylationProfiles,narrowPeakProfiles,pileupCoverageProfiles, progress,last,fac2,genomePath,model,startStrand,endStrand,strandPenalty,nd,t,kmer,use,starts);
				progress.setCurrent(last + 1.0*fac2);
				last += 1.0*fac2;
			}

			ListResult lr = toListResult(ll,talName);

			talRess.add(lr);
			
			protocol.append("...finished predicting "+ll.size()+" sites.\n\n");

		}
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(talRess), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	private String[][] readTALs(FileRepresentation fileContents) throws IOException {
		BufferedReader read = new BufferedReader(new StringReader(fileContents.getContent()));
		
		LinkedList<String[]> tals = new LinkedList<>();
		
		String[] temp = null;
		String str = null;
		while( (str = read.readLine()) != null ){
			if(str.startsWith(">")){
				temp = new String[2];
				temp[1] = str.substring(1).trim().replaceAll("\\s.*", "");
				temp[0] = "";
				tals.add(temp);
			}else{
				tals.getLast()[0] += str;
			}
		}
		read.close();
		
		return tals.toArray(new String[0][]);
		
	}

	@Override
	public String getToolName() {
		return "PrediTALE";
	}

	@Override
	public String getToolVersion() {
		return "0.2";
	}

	@Override
	public String getShortName() {
		return "preditale";
	}

	@Override
	public String getDescription() {
		return "predicts TALE target boxes using a novel model learned from quantitative data and includes epigenetic features";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( QuickTBSPredictionToolMethylAccessibilityAnnotation_fai.class.getClassLoader().getResourceAsStream( "projects/tals/prediction/PrediTALE.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	private void getSites(LinkedList<ComparableElement<ResultSet,Double>> ll, Sequence rvds, String talName,boolean useMethylationData,MethylationprofilHashMap methylationProfiles,NarrowpeakprofilHashMap narrowPeakProfiles,PileupCoverageprofilHashMap pileupCoverageProfiles, ProgressUpdater progress, double last, double fac, String file, QuickScanningSequenceScore model, int startStrand, int endStrand, double strandPenalty, NormalDist nd, double threshold, int kmer, boolean[][] use, int... offs) throws Exception {
		
		BufferedReader read = new BufferedReader(new FileReader(file));
		StringBuffer lastHeader = new StringBuffer();

		long approxTotal = (new File(file)).length();

		int[] pow = new int[kmer];
		int a = (int)model.getAlphabetContainer().getAlphabetLengthAt(0);
		pow[pow.length-1]=1;
		for( int i = pow.length-2; i >= 0; i-- ) {
			pow[i] = pow[i+1]*a;
		}

		CategoricalResult talRes = new CategoricalResult("TALE", "", talName);
		CategoricalResult rvdsRes = new CategoricalResult("RVDs", "", rvds.toString("-", 0, rvds.getLength()));
		
		int[] idxs = new int[offs.length];

		Pair<IntList,ArrayList<Sequence>> pair = null;

		double prog = 0.3;
		Methylationprofil MP = null;
		
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, model.getLength()) ) != null ){
			IntList starts = pair.getFirstElement();
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			int itIdx = 0;
			
			String id_old="";
			
			while( it.hasNext() ) {
				//System.out.println("id_old: "+id_old);
				Sequence seq = it.next();
				
				int sl = seq.getLength();
				int ml = model.getLength();
				prog += (seq.getLength()/(double)approxTotal)*0.7;
				progress.setCurrent(last + prog*fac);
				
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
				if(!(id.equals(id_old))){
					MP=methylationProfiles.getMethylationprofil(id);
					id_old=id;
				}
				
				seq=seq.annotate(true, new MethylationSequenceAnnotation( "methyl", MP ) );
				
				int off = starts.get(itIdx); 

				itIdx++;

				double sPen = 0.0;

				for(int d=startStrand;d<endStrand;d++){
					
					if(d == 1){
						sPen = strandPenalty;
					}
					Arrays.fill(idxs, 0);
					for(int i=0;i<offs.length;i++){
						for(int j=0;j<kmer-1;j++){
							idxs[i] += pow[j+1]*seq.discreteVal(offs[i]+j);
						}
					}

					for(int j=0;j<sl-ml+1;j++){
						//for(int j=0;j<Math.min(sl-ml+1,6);j++){
						for(int i=0;i<idxs.length;i++){
							idxs[i] = (idxs[i]%pow[0])*4 + seq.discreteVal(offs[i]+j+kmer-1);
						}

						boolean used = true;
						for(int k=0;used && k<use.length;k++){
							used &= use[k][idxs[k]];
						}
						if(used){
							
							if(d == 1){
								MP.setStrand(false);
							}else{
								MP.setStrand(true);
							}

							MP.setStartPos(off+j);	

							double score = model.getLogScoreFor(seq, j) + sPen;

							if(score > threshold){
								 ArrayList<Result> rl=new ArrayList<Result>();
								 
								rl.add(new CategoricalResult("Seq-ID", "", id));
								rl.add(new NumericalResult("Position", "",(d==0 ? off+j : off+sl-j-ml) ));
								rl.add(new CategoricalResult("Strand","",d==0 ? "+" : "-"));
								rl.add(new NumericalResult("Score", "", score));
								rl.add(new CategoricalResult("Sequence", "", seq.toString(j, j+model.getLength())));
								rl.add(new NumericalResult("Approx. p-value", "", (1.0-nd.cdf(score))));
								rl.add(rvdsRes);
								rl.add(talRes);
								if(useMethylationData){
									rl.add(new CategoricalResult("MethylationProp", "", ((PFMWrapperTrainSMMethyl) model).getMethylProb(seq,j,j+rvds.getLength())));
								}
								if(narrowPeakProfiles!=null){
									rl.add(new CategoricalResult("isPeakSurroundBox","",narrowPeakProfiles.getNarrowpeakprofil(id).isPeakSurroundPos((d==0 ? off+j : off+sl-j-ml), (d==0 ? true : false))));
								}
								if(pileupCoverageProfiles!=null){
									rl.add(new NumericalResult("countCovPos","",pileupCoverageProfiles.getPileupCoverageprofil(id).getnumberOfCoveragePositionsSurroundPos((d==0 ? off+j : off+sl-j-ml), (d==0 ? true : false))));
								}
									
								 Result[] ra=new Result[rl.size()];
								 rl.toArray(ra);
								ResultSet rs = new ResultSet(ra);
								
								ll.add(new ComparableElement<ResultSet, Double>(rs, -score));
								
								
								
							}
						}

					}
					seq = seq.reverseComplement();
					if(methylationProfiles!=null){
						seq=seq.annotate(true, new MethylationSequenceAnnotation( "methyl", MP ) );
					}
				}
			}
		}
		
		//System.out.println("getSitesEnde");
	}
	
	private TextResult toTextResult(LinkedList<ComparableElement<StringBuffer, Double>>ll, String talName){
		ComparableElement<StringBuffer, Double>[] rsa = ll.toArray(new ComparableElement[0]);
		Arrays.sort(rsa);
		
		StringBuffer sb = new StringBuffer();
		
		sb.append("#Seq-ID\tPosition\tStrand\tScore\tSequence\tApprox. p-value\tRVDs\tTALE\tMethylProb\n");
		for(int i=0;i<rsa.length;i++){
			sb.append(rsa[i].getElement());
			sb.append("\n");
		}
		
		return new TextResult("Predicted binding sites for "+talName, "Predicted binding sites", new FileRepresentation("", sb.toString()), "tsv", this.getClass().getSimpleName(), null, true);
	}
	
	private ListResult toListResult(LinkedList<ComparableElement<ResultSet,Double>> ll, String talName){
		ComparableElement<ResultSet, Double>[] rsa = ll.toArray(new ComparableElement[0]);
		Arrays.sort(rsa);

		ResultSet[] rss = new ResultSet[rsa.length];
		for(int i=0;i<rss.length;i++){
			rss[i] = rsa[i].getElement();
		}
		ListResult lr = new ListResult("Predicted binding sites for "+talName, "Predicted binding sites", null,rss );
		lr.setExport(true);
		return lr;
		
	}

	private NormalDist getThreshold(String file, QuickScanningSequenceScore model2, double p, int startStrand, int endStrand, double strandPenalty) throws Exception{
		BufferedReader read = new BufferedReader(new FileReader(file));
		StringBuffer lastHeader = new StringBuffer();

		Random r = new Random(113);

		p /= 2.0;

		DoubleList scores = new DoubleList();

		Pair<IntList,ArrayList<Sequence>> pair = null;

		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, model2.getLength()) )!= null ){
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			while( it.hasNext() ) {
				Sequence seq = it.next();

				double sPen = 0.0;
				
				for(int d=startStrand;d<endStrand;d++){

					if(d == 1){
						sPen = strandPenalty;
					}
					
					if( p < 1.0){
						double num = p*(seq.getLength()-model2.getLength()+1);
						double meanStep = (seq.getLength()-model2.getLength()+1)/num;
						double sd = Math.sqrt(meanStep);

						for(int j=Math.max(1, (int)Math.round( meanStep + r.nextGaussian()*sd ));j<seq.getLength()-model2.getLength()+1;j+= Math.max(1, (int)Math.round( meanStep + r.nextGaussian()*sd )) ){
							if(j>seq.getLength()-model2.getLength()+1){
								break;
							}
							double score = model2.getLogScoreFor(seq, j) + sPen;

							scores.add(score);
						}
					}else{
						for(int j=0;j<seq.getLength()-model2.getLength()+1;j++ ){

							double score = model2.getLogScoreFor(seq, j) + sPen;

							scores.add(score);
						}
					}
					seq = seq.reverseComplement();
				}
			}
		}
		read.close();
		
		double mean = scores.mean(0, scores.length());
		double meansq = 0;
		double n = 0;

		for(int i=0;i<scores.length();i++){
			double score = scores.get(i); 
			if(score>=mean){
				meansq += score*score;
				meansq += (2*mean-score)*(2*mean-score);
				n += 2;
			}
		}
		meansq /=n;

		double sd = Math.sqrt( meansq - mean*mean );

		NormalDist nd = new NormalDist(mean, sd);

		return nd; 
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String[] getReferences() {
		// TODO Auto-generated method stub
		return null;
	}


}
