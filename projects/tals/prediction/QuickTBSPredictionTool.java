package projects.tals.prediction;

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
import de.jstacs.sequenceScores.statisticalModels.trainable.PFMWrapperTrainSM;
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
import projects.tals.linear.LFModularConditional9C;
import umontreal.ssj.probdist.NormalDist;

public class QuickTBSPredictionTool implements JstacsTool {

	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new QuickTBSPredictionTool());
		cli.run(args);
	}


	public QuickTBSPredictionTool() {

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
			
			
			return new ToolParameterSet(getShortName(),genome,bg,tsel,tals,strand);

		} catch (ParameterException e) {
			throw new RuntimeEOFException(e);
		}


	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {

		progress.setLast(1.0);
		progress.setCurrent(0.0);

		//String modelPath = parameters.getParameterAt(0).getValue().toString();
		boolean backgroundSet = ((SelectionParameter)parameters.getParameterAt(1)).getSelected()==1;
		String backgroundPath = backgroundSet ? ( (ParameterSet)parameters.getParameterAt(1).getValue()).getParameterAt(0).getValue().toString() : parameters.getParameterAt(0).getValue().toString();
		String genomePath = parameters.getParameterAt(0).getValue().toString();
		
		int startStrand = 0, endStrand = 2;
		double strandPenaltyPerc = 0.0;
		if(((SelectionParameter)parameters.getParameterAt(4)).getSelected() == 1){
			endStrand = 1;
		}else if(((SelectionParameter)parameters.getParameterAt(4)).getSelected() == 2){
			startStrand = 1;
		}else{
			strandPenaltyPerc = (double) ((ParameterSet)parameters.getParameterAt(4).getValue()).getParameterAt(0).getValue();
		}
		
		String[][] tals = readTALs(((FileParameter)parameters.getParameterAt(3)).getFileContents());

		LinkedList<Result> talRess = new LinkedList<>();
		
		double fac = 1.0/tals.length;
		double last = 0.0;
		
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
				//System.out.println(currS);
				RVDSequence eff = new RVDSequence(alphabet12, alphabet13, currS);



				//System.out.println("subsamp: "+subsamp);

				//LFModularConditional9C lfmod = new LFModularConditional9C(FileManager.readFile(modelPath));
				LFModularConditional9C lfmod = new LFModularConditional9C( FileManager.readInputStream( QuickTBSPredictionTool.class.getClassLoader().getResourceAsStream( "projects/tals/prediction/preditale_quantitative_PBM.xml" ) ) );
				
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

				//	RVDSequence eff = new RVDSequence(alphabet12, alphabet13, rvdStr);

				QuickScanningSequenceScore model = new PFMWrapperTrainSM(DNAAlphabetContainer.SINGLETON, "", pwm);

				int kmer = Math.min(10, model.getLength()*2/3);
				int[] starts = new int[]{0,(model.getLength()-kmer)/3,(model.getLength()-kmer)*2/3,model.getLength()-kmer};
				protocol.append("Target site length: "+model.getLength()+"\n");
				protocol.append("Using "+kmer+"-mers starting at positions "+Arrays.toString(starts)+"\n");

				NormalDist nd = getThreshold(backgroundPath,model,subsamp,startStrand,endStrand,strandPenalty); 
				progress.setCurrent(last + 0.3*fac2);

				double t = nd.inverseF(1.0-p_value);

				protocol.append("Effective threshold: "+t+"\n");



				boolean[][] use = model.getInfixFilter(kmer, t, starts);


				int n = 0;
				double[] us = new double[use.length];

				for(int i=0;i<use.length;i++){
					for(int j=0;j<use[i].length;j++,n++){
						if(use[i][j]){
							us[i]++;
						}
					}

				}
				protocol.append("Number of "+kmer+"-mers passing filtering: "+Arrays.toString(us)+"\n\n");
				//System.out.println(Arrays.toString(us));
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
				getSites(ll, eff, talName, progress,last,fac2,genomePath,model,startStrand,endStrand,strandPenalty,nd,t,kmer,use,starts);
				progress.setCurrent(last + 1.0*fac2);
				last += 1.0*fac2;

			}

			ListResult lr = toListResult(ll,talName);
			//TextResult lr = toTextResult(ll, talName);

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
				if(tals.getLast()[0].indexOf("-")<0) {
					throw new IOException("Malformed TALE sequence; use dashes to separate RVDs");
				}
				if(tals.getLast()[0].length()>500) {
					throw new IOException("Malformed TALE sequence; more than 150 RVDs per TALE not allowed");
				}
				
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
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "preditale";
	}

	@Override
	public String getDescription() {
		return "predicts TALE target boxes using a novel model learned from quantitative data";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( QuickTBSPredictionTool.class.getClassLoader().getResourceAsStream( "projects/tals/prediction/PrediTALE.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}







	private void getSites(LinkedList<ComparableElement<ResultSet,Double>> ll, Sequence rvds, String talName, ProgressUpdater progress, double last, double fac, String file, QuickScanningSequenceScore model, int startStrand, int endStrand, double strandPenalty, NormalDist nd, double threshold, int kmer, boolean[][] use, int... offs) throws Exception {
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

		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, model.getLength()) ) != null ){
			IntList starts = pair.getFirstElement();
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			int itIdx = 0;
			
			

			while( it.hasNext() ) {
				Sequence seq = it.next();

				int sl = seq.getLength();
				int ml = model.getLength();
				String rvdStr = rvds.toString("-", 0, rvds.getLength());
				
				prog += (seq.getLength()/(double)approxTotal)*0.7;
				progress.setCurrent(last + prog*fac);

				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
				int off = starts.get(itIdx);
				itIdx++;

				double sPen = 0.0;
				
				for(int d=startStrand;d<endStrand;d++){//TODO revcom if startStrand==1

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
						for(int i=0;i<idxs.length;i++){
							idxs[i] = (idxs[i]%pow[0])*4 + seq.discreteVal(offs[i]+j+kmer-1);
						}


						boolean used = true;
						for(int k=0;used && k<use.length;k++){
							used &= use[k][idxs[k]];
						}
						if(used){
							double score = model.getLogScoreFor(seq, j) + sPen;

							if(score > threshold){
								
								/*StringBuffer sb = new StringBuffer();
								sb.append(id);
								sb.append("\t");
								sb.append((d==0 ? off+j : off+sl-j-ml));
								sb.append("\t");
								sb.append(d==0 ? "+" : "-");
								sb.append("\t");
								sb.append(score);
								sb.append("\t");
								sb.append(seq.toString(j, j+model.getLength()));
								sb.append("\t");
								sb.append((1.0-nd.cdf(score)));
								sb.append("\t");
								sb.append(rvdStr);
								sb.append("\t");
								sb.append(talName);
								
								ll.add(new ComparableElement<StringBuffer, Double>(sb, -score));*/
								
								ResultSet rs = new ResultSet(new Result[]{
										new CategoricalResult("Seq-ID", "", id),
										new NumericalResult("Position", "",(d==0 ? off+j : off+sl-j-ml) ),
										new NumericalResult("Distance to end", "",seq.getLength()-(d==0 ? off+j+ml : off+sl-j) ),
										new CategoricalResult("Strand","",d==0 ? "+" : "-"),
										new NumericalResult("Score", "", score),
										new CategoricalResult("Sequence", "", seq.toString(j, j+model.getLength())),
										new NumericalResult("Approx. p-value", "", (1.0-nd.cdf(score))),
										rvdsRes,
										talRes
										
								});
								ll.add(new ComparableElement<ResultSet, Double>(rs, -score));
								
								
								
							}
						}

					}



					seq = seq.reverseComplement();
				}

			}
			

		}
		

	}
	
	private TextResult toTextResult(LinkedList<ComparableElement<StringBuffer, Double>>ll, String talName){
		ComparableElement<StringBuffer, Double>[] rsa = ll.toArray(new ComparableElement[0]);
		Arrays.sort(rsa);
		
		StringBuffer sb = new StringBuffer();
		
		sb.append("#Seq-ID\tPosition\tStrand\tScore\tSequence\tApprox. p-value\tRVDs\tTALE\n");
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

		//	System.out.println("meansq: "+meansq+" n: "+n);

		meansq /=n;

		double sd = Math.sqrt( meansq - mean*mean );


		NormalDist nd = new NormalDist(mean, sd);
			//System.out.println("mean: "+mean+" sd: "+sd);
		return nd; 
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
		return new String[]{
				"@article{erkes19preditale,\n" + 
				"	title = {{PrediTALE}: A novel model learned from quantitative data allows for new perspectives on {TALE} targeting},\n" + 
				"	author = {Erkes, Annett AND M\\\"ucke, Stefanie AND Reschke, Maik AND Boch, Jens AND Grau, Jan},\n" + 
				"	journal = {PLOS Computational Biology},\n" + 
				"	year = {2019},\n" + 
				"	volume = {15},\n" + 
				"	number = {7},\n" + 
				"	pages = {1-31},\n" + 
				"	doi = {10.1371/journal.pcbi.1007206}\n" + 
				"	}\n"
		};
	}


}
