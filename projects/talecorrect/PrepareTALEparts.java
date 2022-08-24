package projects.talecorrect;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.data.DNADataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
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

public class PrepareTALEparts implements JstacsTool{

	public static void main(String[] args) throws Exception {
		CLI cli = new CLI(new PrepareTALEparts());
		cli.run(args);
	}

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();

		try {
			pars.add(new FileParameter("TALE DNA parts", "The TALE DNA parts output file of AnnoTALE analyze.", "fasta,fa,fas,fna", true));
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
		return new ToolParameterSet(this.getShortName(), pars.toArray(new Parameter[0]));	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol,
			ProgressUpdater progress, int threads) throws Exception {
		FileRepresentation TALEdnaPartsFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents();
		
		SimpleSequenceAnnotationParser parser = new SimpleSequenceAnnotationParser();
		
		DNADataSet ds = new DNADataSet( TALEdnaPartsFile.getFilename(),'>', parser );

		StringBuffer NSB=new StringBuffer();
		StringBuffer RSB=new StringBuffer();
		StringBuffer CSB=new StringBuffer();
		String header="";
		String nextheader="";
		for(int i=0;i<ds.getNumberOfElements();i++){
			header=ds.getElementAt(i).getSequenceAnnotationByType("unparsed comment line", 0).getResultAt(0).getValue().toString();
			//System.out.println(header);
			if(header.matches(".*N-terminus.*")){
				nextheader=ds.getElementAt(i+1).getSequenceAnnotationByType("unparsed comment line", 0).getResultAt(0).getValue().toString();
				if(nextheader.matches(".*repeat 1.*")){
					NSB.append(">"+header+"\n");
					NSB.append(ds.getElementAt(i).toString()+ds.getElementAt(i+1).toString(0, 10)+"\n");
				}
				
			}
			if(header.matches(".*repeat.*")){
				RSB.append(">"+header+"\n");
				RSB.append(ds.getElementAt(i).toString()+"\n");
			}
			if(header.matches(".*C-terminus.*")){
				CSB.append(">"+header+"\n");
				CSB.append(ds.getElementAt(i).toString()+"\n");
			}
		}	

		protocol.append( "\nWriting outputs.\n" );
		
		TextResult frN = new TextResult( "TALE_DNA_parts.N-terminus.10bpRepeat1", "Output with N-terminus sequences including 10 bp of repeat 1.", new FileRepresentation( "", NSB.toString() ), "fasta", getToolName(), "fasta/dna",true );
		TextResult frR = new TextResult( "TALE_DNA_parts.repeat", "Output with repeat sequences.", new FileRepresentation( "", RSB.toString() ), "fasta", getToolName(), "fasta/dna",true );
		TextResult frC = new TextResult( "TALE_DNA_parts.C-terminus", "Output with C-terminus sequences.", new FileRepresentation( "", CSB.toString() ), "fasta", getToolName(), "fasta/dna",true );
		
		ResultSet set = new ResultSet( new Result[]{frN,frR,frC} );
		return new ToolResult("Result of "+getToolName(), getToolName()+" on \""+TALEdnaPartsFile.getFilename()+"\"", null, set, parameters, getToolName(), new Date(System.currentTimeMillis()) );

	}

	@Override
	public String getToolName() {
		return "PrepareTALEparts";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "prepare";
	}

	@Override
	public String getDescription() {
		// TODO Auto-generated method stub
		return "Prepares output of AnnoTALE alayze - TALE dna parts file for hmmbuild.";
	}

	@Override
	public String getHelpText() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		// TODO Auto-generated method stub
		return null;
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
