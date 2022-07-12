package projects.talecorrect;

import java.util.List;
import java.util.Collection;
import java.util.Collections;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedList;
import java.io.BufferedReader;
import java.io.StringReader;
import de.jstacs.data.sequences.Sequence;
import java.util.HashMap;
import de.jstacs.data.DataSet;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;

public class PolishTALESubstitutions implements JstacsTool {
	public static void main( String[] args) throws Exception {
		CLI cli = new CLI(new PolishTALESubstitutions());
		cli.run(args);
	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();

		try {
			pars.add(new FileParameter("Corrected Assembly", "The corrected assembly - output file of CorrectTALESequences.", "fasta,fa,fas,fna", true));
			pars.add(new FileParameter("IGVtoolsCountWig", "The output of igvtools count - wig-File.", "wig", true));
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
		return new ToolParameterSet(this.getShortName(), pars.toArray(new Parameter[0]));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		FileRepresentation nanoAssemblyTALECorrectedFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents();
		FileRepresentation inputIGVtoolsCountWig = ((FileParameter)parameters.getParameterAt(1)).getFileContents();
		
		SimpleSequenceAnnotationParser parser = new SimpleSequenceAnnotationParser();
		String[] symbols_caseSensitive = { "A", "C", "G", "T", "a", "c", "g", "t" };
		 AlphabetContainer conCaseSensitive = new AlphabetContainer(new DiscreteAlphabet(false, symbols_caseSensitive));
		 
		 DataSet ds = new DataSet(conCaseSensitive, new SparseStringExtractor(nanoAssemblyTALECorrectedFile.getFilename(), '>',parser));
		 AlphabetContainer con = ds.getAlphabetContainer();
		 HashMap<String, Sequence> seqHM = new HashMap<String, Sequence>();
		for (int i = 0; i < ds.getNumberOfElements(); ++i) {
			 Sequence seq = ds.getElementAt(i);
			 String header = seq.getAnnotation()[0].getResultAt(0).getValue().toString();
			 String seqName = header.split(" ")[0];
			seqHM.put(seqName, seq);
		}

		 BufferedReader BR = new BufferedReader(new StringReader(inputIGVtoolsCountWig.getContent()));
		Sequence seqCorrected = null;
		String line = "";
		String chrom = "";
		String[] igvSplitLine = new String[8];
		 ArrayList<Integer> igvCounts = new ArrayList<Integer>();
		while ((line = BR.readLine()) != null) {
			if (line.matches("variableStep chrom.*")) {
				chrom = line.split(" ")[1].substring(6);
			}
			else {
				if (line.matches("track type=.*") || line.matches("#Columns: Pos.*")) {
					continue;
				}
				igvSplitLine = line.split("\t");
				 int assemblyPos = Integer.parseInt(igvSplitLine[0]);
				igvCounts.clear();
				for (int j = 1; j <= 4; ++j) {
					igvCounts.add((int)Double.parseDouble(igvSplitLine[j]));
				}
				 Integer maxVal = Collections.max((Collection<? extends Integer>)igvCounts);
				 Integer maxIdx = igvCounts.indexOf(maxVal);
				 String nuclUnpolished = seqHM.get(chrom).getSubSequence(assemblyPos - 1, 1).toString();
				 String nuclMaxCounts = con.getSymbol(0, maxIdx);
				if (nuclUnpolished.toUpperCase().equals(nuclMaxCounts)) {
					continue;
				}
				int nuclUnpolished_count = igvCounts.get(seqHM.get(chrom).discreteVal(assemblyPos - 1) % 4);
				 double frac = maxVal / 1.5;
			   
				if (frac < nuclUnpolished_count) {
					continue;
				}
				System.out.println("chrom: " + chrom + "\t" + "pos: " + assemblyPos);
				System.out.println(line);
				//System.out.println("chrom: " + chrom + "\t" + "pos: " + assemblyPos);
				System.out.println("nuclUnpolished: " + seqHM.get(chrom).getSubSequence(assemblyPos - 1, 1) + " vs nuclMaxCounts: " + con.getSymbol(0, maxIdx) + "(" + maxVal + ")");
				//System.out.println("(" + maxVal + "/" + 1.5 + "=" + frac + ")>=" + nuclUnpolished_count);
				System.out.println();
				 String seqBefore = seqHM.get(chrom).getSubSequence(0, assemblyPos - 1).toString();
				 String seqAfter = seqHM.get(chrom).getSubSequence(assemblyPos).toString();
				 String substitution = con.getSymbol(0, maxIdx).toLowerCase();
				seqCorrected = Sequence.create(conCaseSensitive, seqBefore + substitution + seqAfter);
				seqHM.replace(chrom, seqCorrected);
			}
		}
		BR.close();
		
		StringBuffer sb = new StringBuffer();
		
		 List<String> sortedKeys = new ArrayList<String>(seqHM.keySet());
		Collections.sort(sortedKeys);
		for ( String seqName : sortedKeys) {
			//System.out.println(seqName);
			sb.append(">" + seqName + "\n");
			 String seq2 = seqHM.get(seqName).toString();
			 sb.append(seq2 + "\n");
		}
		LinkedList<Result> ress = new LinkedList<>();
		TextResult tr = new TextResult( "polishedSubstitutions", "", new FileRepresentation( "", sb.toString() ), "fa", getToolName(), "fasta/as",true );
		ress.addFirst(tr);
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(ress.toArray(new Result[0])), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "PolishTALESubstitutions";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "polish";
	}

	@Override
	public String getDescription() {
		return "polishes substitions after run of CorrectTALESequences tool";
	}

	@Override
	public String getHelpText() {
		// TODO Auto-generated method stub
		return null;
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
		// TODO Auto-generated method stub
		return null;
	}
}