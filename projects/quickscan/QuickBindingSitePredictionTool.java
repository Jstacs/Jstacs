package projects.quickscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.zip.GZIPInputStream;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
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
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.LargeSequenceReader;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.dimont.ThresholdedStrandChIPper;
import umontreal.iro.lecuyer.probdist.NormalDist;

public class QuickBindingSitePredictionTool implements JstacsTool {

	public static void main(String[] args) throws Exception{
		CLI cli = new CLI(new QuickBindingSitePredictionTool());
		cli.run(args);
	}
	
	
	public QuickBindingSitePredictionTool() {
	
	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter model = new FileParameter("Dimont model", "The model returned by Dimont (in XML format)", "xml", true);
		
		FileParameter genome = new FileParameter("Sequences","The sequences (e.g., a genome) to scan for binding sites","fa,fas,fasta,fa.gz,fas.gz,fasta.gz",true);
		
		
		FileParameter background = new FileParameter("Background sequences","The sequences (e.g., a genome) for determining the prediction threshold","fa,fas,fasta,fa.gz,fas.gz,fasta.gz",true);
		
		
		
		try {
			 
			SelectionParameter bg = new SelectionParameter(DataType.PARAMETERSET,new String[]{"sub-sample","background sequences"}, new ParameterSet[]{
					new SimpleParameterSet(),
					new SimpleParameterSet(background)
			},"Background sample","The sequences for determining the prediction threshold. Either a sub-sample of the input sequences or a dedicated background data set.",true);
			
			SimpleParameter pvalue = new SimpleParameter(DataType.DOUBLE, "Significance level", "The significance level for determining the prediction threshold", true, new NumberValidator<Comparable<Double>>(0.0, 1E-3), 1E-6);
			
			SimpleParameter num = new SimpleParameter(DataType.INT,"Number of sites", "The number of expected binding sites for determining the prediction threshold", true, new NumberValidator<Comparable<Integer>>(1, 1000000),10000);
			
			SelectionParameter tsel = new SelectionParameter(DataType.PARAMETERSET, new String[]{"significance level","number of sites"}, new ParameterSet[]{
					new SimpleParameterSet(pvalue),
					new SimpleParameterSet(num)
			}, "Threshold specification", "The way of defining the prediction threshold. Either by explicitly defining a significance level or by specifying the number of expected sites", true);
			
			
			return new ToolParameterSet(getShortName(),model,genome,bg,tsel);
			
		} catch (ParameterException e) {
			throw new RuntimeException(e);
		}
		
		
	}
	
	private static long getUncompressedSize(FileRepresentation genomeFile) throws FileNotFoundException, IOException	{
		if((new File(genomeFile.getFilename()).exists())){
			String genomePath = genomeFile.getFilename();
			if(genomePath.endsWith(".gz")) {
				long size = -1;
				try (RandomAccessFile fp = new RandomAccessFile(new File(genomePath), "r")) {        
					fp.seek(fp.length() - Integer.BYTES);
					int n = fp.readInt();
					size = Integer.toUnsignedLong(Integer.reverseBytes(n));
				}
				return size;
			}else {
				return (new File(genomePath)).length();
			}
		}else {
			return genomeFile.getContent().length();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		progress.setLast(1.0);
		progress.setCurrent(0.0);
		
		FileRepresentation modelFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents();
		boolean backgroundSet = ((SelectionParameter)parameters.getParameterAt(2)).getSelected()==1;
		FileRepresentation backgroundFile = ((FileParameter)(backgroundSet ? ( (ParameterSet)parameters.getParameterAt(2).getValue()).getParameterAt(0) : parameters.getParameterAt(1))).getFileContents();
		FileRepresentation genomeFile = ((FileParameter) parameters.getParameterAt(1)).getFileContents();
		
		protocol.appendHeading("Starting predictions...\n");
		protocol.append("Using "+(backgroundSet? " background set.\n" : " sub-sample of input data.\n"));
		
		boolean byp = ((SelectionParameter)parameters.getParameterAt(3)).getSelected()==0;
		double p_value = 0;
		if(byp){
			p_value = (Double) ((ParameterSet)parameters.getParameterAt(3).getValue()).getParameterAt(0).getValue();
		}else{
			int nsites = (Integer) ((ParameterSet)parameters.getParameterAt(3).getValue()).getParameterAt(0).getValue();
			p_value = nsites/(double)getUncompressedSize(genomeFile)/2.0;
		}
		protocol.append("Significance level: "+p_value+"\n");
		protocol.appendWarning("The p-values and, hence, the significance level are only approximate values and may not fully reflect the number of predictions for a specific input file.\n");
		
		double subsamp = 0;
		if(backgroundSet){
			subsamp = 1;
		}else{
			subsamp = 1E6/(double)getUncompressedSize(backgroundFile);
		}
		//System.out.println("subsamp: "+subsamp);
		
		GenDisMixClassifier cl = new GenDisMixClassifier(new StringBuffer(modelFile.getContent()));
		
		ThresholdedStrandChIPper fg =(ThresholdedStrandChIPper) cl.getDifferentiableSequenceScore(0);
		QuickScanningSequenceScore lslim = (QuickScanningSequenceScore) fg.getFunction(0);
		
		int kmer = 10;
		if(kmer > lslim.getLength()) {
			kmer = (lslim.getLength()+1)/2;
		}
		int[] starts = new int[]{0,(lslim.getLength()-kmer)/2,lslim.getLength()-kmer};
		
		NormalDist nd = getThreshold(backgroundFile,lslim,subsamp); 
		progress.setCurrent(0.3);
		
		double t = nd.inverseF(1.0-p_value);

		protocol.append("Effective threshold: "+t+"\n");
		
		
		
		boolean[][] use = lslim.getInfixFilter(kmer, t, starts);
		
		
		int n = 0;
		double[] us = new double[use.length];
		
		for(int i=0;i<use.length;i++){
			for(int j=0;j<use[i].length;j++,n++){
				if(use[i][j]){
					us[i]++;
				}
			}
			
		}
		
		int[] o = ToolBox.order(us, false);
		boolean[][] temp = new boolean[use.length][];
		int[] temps = new int[starts.length];
		for(int i=0;i<use.length;i++){
			temp[i] = use[o[i]];
			temps[i] = starts[o[i]];
		}
		use = temp;
		starts = temps;
		
		protocol.appendHeading("Predicting sites...\n");
		Pair<Integer,TextResult> tr = getSites(progress,genomeFile,lslim,nd,t,kmer,use,starts);
		progress.setCurrent(1.0);
		protocol.append("...finished predicting "+tr.getFirstElement()+" sites.\n");
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr.getSecondElement()), parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}

	@Override
	public String getToolName() {
		return "Quick Prediction Tool";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "quickpred";
	}

	@Override
	public String getDescription() {
		return "predicts binding sites for a fixed threshold";
	}

	@Override
	public String getHelpText() {
		return "**Quick Prediction Tool** predicts binding sites of a transcription factor based on a motif model and is also suited for genome-wide predictions. The motif model is provided as the XML output of (Slim) Dimont. "
				+ "\n\n"
				+ "The tool outputs a list of predictions including, for every prediction, the ID"
				+ "of the sequence (e.g., chromosome) containing the binding site, position and strand of the matching sub-sequence, its score according to the model, the"
				+ " sub-sequence itself (in strand orientation according to the model), and a p-value from a normal distribution fitted to the score distribution of the provided negative examples "
				+ "or a sub-sample of the input data (parameter \"Background sample\")."
				+ "\n\n"
				+ "If you experience problems using Quick Prediction Tool, please contact_ us.\n"
				+ ".. _contact: mailto:grau@informatik.uni-halle.de";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	
	
	
	
	
	
	private Pair<Integer,TextResult> getSites(ProgressUpdater progress, FileRepresentation fileRep, QuickScanningSequenceScore model, NormalDist nd, double threshold, int kmer, boolean[][] use, int... offs) throws Exception {
		BufferedReader read = null;
		if((new File(fileRep.getFilename())).exists()) {
			String file = fileRep.getFilename();
			if(file.toLowerCase().endsWith(".gz")) {
				read = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			}else {
				read = new BufferedReader(new FileReader(file));
			}
		}else {
			read = new BufferedReader(new StringReader(fileRep.getContent()));
		}
		StringBuffer lastHeader = new StringBuffer();
		
		long approxTotal = getUncompressedSize(fileRep);
		
		int[] pow = new int[kmer];
		int a = (int)model.getAlphabetContainer().getAlphabetLengthAt(0);
		pow[pow.length-1]=1;
		for( int i = pow.length-2; i >= 0; i-- ) {
			pow[i] = pow[i+1]*a;
		}
		
		int nCorr = 0;
		int nWro = 0;
		long nTot = 0;
		
		
		int[] idxs = new int[offs.length];
		
		StringBuffer sb = new StringBuffer();
		sb.append("Seq-ID\tPosition\tStrand\tScore\tSequence\tApprox. p-value\n");
		
		Pair<IntList,ArrayList<Sequence>> pair = null;
		
		double prog = 0.3;
		
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, model.getLength(), model.getAlphabetContainer()) ) != null ){
			IntList starts = pair.getFirstElement();
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			int itIdx = 0;
			
			
			while( it.hasNext() ) {
				Sequence seq = it.next();
				
				prog += (seq.getLength()/(double)approxTotal)*0.7;
				progress.setCurrent(prog);
				
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
				int off = starts.get(itIdx);
				itIdx++;
				
				for(int d=0;d<2;d++){

					Arrays.fill(idxs, 0);
					for(int i=0;i<offs.length;i++){
						for(int j=0;j<kmer-1;j++){
							idxs[i] += pow[j+1]*seq.discreteVal(offs[i]+j);
						}
					}
					
					
					for(int j=0;j<seq.getLength()-model.getLength()+1;j++,nTot++){
						for(int i=0;i<idxs.length;i++){
							idxs[i] = (idxs[i]%pow[0])*a + seq.discreteVal(offs[i]+j+kmer-1);
						}
						
						
						boolean used = true;
						for(int k=0;used && k<use.length;k++){
							used &= use[k][idxs[k]];
						}
						if(used){
							double score = model.getLogScoreFor(seq, j);
							
							if(score > threshold){
								
								sb.append(id+"\t"+(d==0 ? off+j : off+seq.getLength()-j-model.getLength())+"\t"+(d==0 ? "+" : "-")+"\t"+score+"\t"+seq.toString(j, j+model.getLength())+"\t"+(1.0-nd.cdf(score))+"\n");
								/*
								 * ResultSet rs = new ResultSet(new Result[]{ new CategoricalResult("Seq-ID",
								 * "", id), new NumericalResult("Position", "",(d==0 ? off+j :
								 * off+seq.getLength()-j-model.getLength()) ), new
								 * CategoricalResult("Strand","",d==0 ? "+" : "-"), new NumericalResult("Score",
								 * "", score), new CategoricalResult("Sequence", "", seq.toString(j,
								 * j+model.getLength())), new NumericalResult("Approx. p-value", "",
								 * (1.0-nd.cdf(score))) }); 
								 * ll.add(rs);
								 */
								nCorr++;
							}else{
								nWro++;
							}
						}
						
					}
					
					
					
					seq = seq.reverseComplement();
				}
				
			}
			
			
		}
	//	System.out.println("nTot: "+nTot);
		//return new ListResult("Predicted binding sites", "Predicted binding sites for threshold "+threshold, null,ll );
		return new Pair<Integer,TextResult>(nCorr,new TextResult("Predicted binding sites", "Predicted binding sites for threshold "+threshold, new FileParameter.FileRepresentation("", sb.toString()),"tsv","QBSPT",null,true));
		
	}
	
	private NormalDist getThreshold(FileRepresentation fileRep, QuickScanningSequenceScore model2, double p) throws Exception{
		BufferedReader read = null;
		if((new File(fileRep.getFilename())).exists()) {
			String file = fileRep.getFilename();
			if(file.toLowerCase().endsWith(".gz")) {
				read = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			}else {
				read = new BufferedReader(new FileReader(file));
			}
		}else {
			read = new BufferedReader(new StringReader(fileRep.getContent()));
		}
		StringBuffer lastHeader = new StringBuffer();
		
		
		Random r = new Random(113);
		
		p /= 2.0;
		
		DoubleList scores = new DoubleList();
		
		Pair<IntList,ArrayList<Sequence>> pair = null;
		
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, model2.getLength(), model2.getAlphabetContainer()) )!= null ){
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			while( it.hasNext() ) {
				Sequence seq = it.next();
				
				for(int d=0;d<2;d++){

					if( p < 1.0){
						double num = p*(seq.getLength()-model2.getLength()+1);
						double meanStep = (seq.getLength()-model2.getLength()+1)/num;
						double sd = Math.sqrt(meanStep);

						for(int j=Math.max(1, (int)Math.round( meanStep + r.nextGaussian()*sd ));j<seq.getLength()-model2.getLength()+1;j+= Math.max(1, (int)Math.round( meanStep + r.nextGaussian()*sd )) ){
							if(j>seq.getLength()-model2.getLength()+1){
								break;
							}
							double score = model2.getLogScoreFor(seq, j);

							scores.add(score);
						}
					}else{
						for(int j=0;j<seq.getLength()-model2.getLength()+1;j++ ){
							
							double score = model2.getLogScoreFor(seq, j);

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
	//	System.out.println("mean: "+mean+" sd: "+sd);
		return nd; 
	}

	
	@Override
	public ToolResult[] getTestCases(String path) {
		try {
			return new ToolResult[]{
					new ToolResult(FileManager.readFile(path+File.separator+"xml/qbspt.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public void clear() {
		
	}

	@Override
	public String[] getReferences() {
		return null;
	}

	
}
