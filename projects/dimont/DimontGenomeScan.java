package projects.dimont;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mixture.StrandDiffSM.InitMethod;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;

/**
 * This class allows to do a genomewide scan using a motif model.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class DimontGenomeScan implements JstacsTool {

	public static void main(String[] args) throws Exception {
		
		CLI cl = new CLI(new DimontGenomeScan());
		
		cl.run(args);
		
	}

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> parameters = new LinkedList<Parameter>();
		
		parameters.add(new FileParameter("Dimont classifier", "The classifier from the Dimont output for one motif", "xml", true));
		
		parameters.add(new FileParameter("Input file", "The file containing the sequences to be scanned (e.g., a genome)", "fasta,fa,fas", true));
		
		try {
			parameters.add(new SimpleParameter(DataType.DOUBLE, "Threshold", "Threshold on the required per-base probability", true, new NumberValidator<Double>(0d, 1d), 0.25 ));
			parameters.add(new SimpleParameter(DataType.BOOLEAN, "Best Strand", "switch which allows to output at a specific position only the best strand or both strands if the corresponding score is above the threshold", true, true ));
		} catch (Exception e) {
			e.printStackTrace();
		}	
		
		return new ToolParameterSet(getShortName(),parameters.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		GenDisMixClassifier cl = new GenDisMixClassifier( new StringBuffer(((FileParameter)parameters.getParameterAt(0)).getFileContents().getContent()) );
		
		DifferentiableStatisticalModel model = ((AbstractSingleMotifChIPper) cl.getDifferentiableSequenceScore(0)).getFunction(0);
		
		StrandDiffSM model2 = new StrandDiffSM(model, 1, true, InitMethod.INIT_FORWARD_STRAND, 0.5);
		
		
		double threshold = Math.log((Double)parameters.getParameterAt(2).getValue())*model.getLength();
		
		StringBuffer lastHeader = new StringBuffer();
		BufferedReader read = new BufferedReader(new FileReader(((FileParameter)parameters.getParameterAt(1)).getFileContents().getFilename()));
		
		File out = File.createTempFile("dimontscan", "_dgs.temp", new File("."));
		out.deleteOnExit(); 
		
		SafeOutputStream sos = SafeOutputStream.getSafeOutputStream(new FileOutputStream(out));

		boolean best = (Boolean) parameters.getParameterAt(3).getValue();
		while( readNextSequences(read, lastHeader, model2.getLength()) ){
			Iterator<Sequence> it = seqs.iterator();
			int i = 0;
			while( it.hasNext() ) {
				Sequence seq = it.next();
				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
				int off = starts.get(i);
				
				for(int j=0;j<seq.getLength()-model2.getLength()+1;j++){
					double[] compScore = model2.getComponentScores(seq, j);
					
					int idx = compScore[0] >= compScore[1] ? 0 : 1;
					
					if( compScore[idx] >= threshold ) {
						if( compScore[0] > threshold && (!best || idx==0) ){
							sos.writeln(id+"\t"+(off+j)+"\t"+compScore[0]+"\t+");
						}
						if( compScore[1] > threshold && (!best || idx==1) ){
							sos.writeln(id+"\t"+(off+j)+"\t"+compScore[1]+"\t-");
						}
					}					
				}				
				i++;				
			}
			
		}
		
		sos.close();
		
		return new ToolResult("Dimont predictions", "", null, new ResultSet( new TextResult("Dimont predictions", "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "txt", getToolName(), null, true)), parameters, getToolName(), new Date());
		
		
	}
	
	
	static IntList starts = new IntList();
	static ArrayList<Sequence> seqs = new ArrayList<Sequence>();
	
	//TODO large chromosomes
	public static boolean readNextSequences(BufferedReader read, StringBuffer lastHeader,  int modelLength ) throws Exception {
		//System.out.println("started reading");
		String str = null;
		
		StringBuffer line = new StringBuffer();
		
		starts.clear();
		seqs.clear();
		
		Pattern acgt = Pattern.compile( "[ACGT]+", Pattern.CASE_INSENSITIVE );
		
		AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
		
		int size = 0;
		
		while( (str = read.readLine()) != null || line.length() > 0 ){
			if(str != null){
				str = str.trim();
			}
			if(str == null || str.startsWith( ">" )){//next sequence
				String header = lastHeader.toString();
				if(str != null){
					lastHeader.delete( 0, lastHeader.length() );
					lastHeader.append( str.substring( 1 ).trim() );
				}
				if(line.length() > 0){//we have a sequence
					int idx = header.indexOf(" ");
					if( idx > 0 ) {//TODO new
						header = header.substring(0, idx);
					}
					SequenceAnnotation annotation = new SequenceAnnotation( "id", header );
					
					String seqStr = line.toString();
					line.delete( 0, line.length() );
					Matcher match = acgt.matcher( seqStr );
					while(match.find()){
						int start = match.start();
						int end = match.end();
						int l = end-start;
						if( l >= modelLength ) {
							Sequence seq = Sequence.create( con, seqStr.substring( start, end ) );
							seq = seq.annotate( false, annotation );
							seqs.add( seq );
							size += l;
							starts.add( start );
							//	ends.add( seqStr.length()-end );
						}
					}
					if(size > 1E7 || str == null){
						return true;
					}
				}
			}else{
				line.append( str );
			}	
		}
		return false;
	}
	
	

	@Override
	public String getToolName() {
		return "Dimont genome scan";
	}

	@Override
	public String getToolVersion() {
		return "1.0";
	}

	@Override
	public String getShortName() {
		return "scan";
	}

	@Override
	public String getDescription() {
		return "scans a genome for prediction of a Dimont model";
	}

	@Override
	public String getHelpText() {
		return "";
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
		return null;
	}

}
