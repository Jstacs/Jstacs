package projects.sigma;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.hmm.models.DifferentiableHigherOrderHMM;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.IntList;
import de.jstacs.utils.LargeSequenceReader;
import de.jstacs.utils.Pair;

public class GenomicScan implements JstacsTool {

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new boolean[]{false}, new GenomicScan());
		cli.run(args);
		
	}
	
	
	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<>();
		
		pars.add(new FileParameter("Input sequences", "", "fasta,fas,fa", true));
		pars.add(new FileParameter("Model", "", "xml", true));
		
		try {
			pars.add(new SimpleParameter(DataType.INT, "length", "Sub-sequence length", true));
		} catch (DatatypeNotValidException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
		
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		FileParameter fp = (FileParameter) parameters.getParameterAt(1);
		FileRepresentation fr = fp.getFileContents();
		
		DifferentiableHigherOrderHMM hmm = new DifferentiableHigherOrderHMM(FileManager.read(new StringReader(fr.getContent())));
		
		String genomePath = parameters.getParameterAt(0).getValue().toString();
		BufferedReader read = new BufferedReader(new FileReader(genomePath));
		
		int length = (int) parameters.getParameterAt(2).getValue();
		
		StringBuffer lastHeader = new StringBuffer();
		
		Pair<IntList,ArrayList<Sequence>> pair = null;
		
		File tempFile = File.createTempFile("scan", ".tsv");
		System.out.println(tempFile);
		PrintWriter wr = new PrintWriter(tempFile);
		
		while( (pair = LargeSequenceReader.readNextSequences(read, lastHeader, length) ) != null ){
			
			IntList starts = pair.getFirstElement();
			ArrayList<Sequence> seqs = pair.getSecondElement();
			Iterator<Sequence> it = seqs.iterator();
			int itIdx = 0;
			
			

			while( it.hasNext() ) {
				Sequence seq = it.next();

				int sl = seq.getLength();
				int ml = length;

				String id = seq.getSequenceAnnotationByType("id", 0).getIdentifier().trim();
				int off = starts.get(itIdx);
			
				for(int d=0;d<2;d++){
					
					for(int j=0;j<sl-ml+1;j++){
						
						Pair<IntList,Double> path = hmm.getViterbiPathFor(j, j+length-1, seq);
						
						wr.println(id+"\t"+(d==0 ? (off+j) : (off+sl-j-ml) )+"\t"+(d==0?"+":"-")+"\t"+path.getSecondElement());
					}
					
					seq = seq.reverseComplement();
				}
				
				
			}
		}
		
		wr.close();
		
		FileRepresentation file = new FileRepresentation(tempFile.getAbsolutePath());
		
		TextResult tr = new TextResult("Profile", "", file, "tsv", this.getShortName(), null, true);
		
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Genomic Scan";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "scan";
	}

	@Override
	public String getDescription() {
		return "";
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
