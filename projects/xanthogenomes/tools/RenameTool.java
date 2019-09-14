package projects.xanthogenomes.tools;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.StringReader;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.biojava.bio.seq.Feature;
import org.biojava.ontology.Term;
import org.biojavax.Note;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequence.IOTools;
import org.biojavax.bio.seq.RichSequenceIterator;

import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.JstacsTool.ResultEntry;


public class RenameTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter names = new FileParameter( "Rename Table", "A tab-separated table containing the old name in the first column and the new name in the second column. Output \"TALE names\" of \"TALE Class Assignment\" tool.", "tsv", true );
		
		FileParameter input = new FileParameter( "Input file", "The input Genbank or GFF3 file that should be renamed.", "gb,gbk,genbank,gff,gff3", true );
		
		return new ToolParameterSet( getShortName(), names, input );
	}

	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		
		FileParameter names = (FileParameter)parameters.getParameterAt( 0 );
		FileParameter input = (FileParameter)parameters.getParameterAt( 1 );
		
		
		HashMap<String, String> nameMap = new HashMap<String,String>();
		
		String[] lines = names.getFileContents().getContent().split( "\n" );
		for(int i=0;i<lines.length;i++){
			lines[i] = lines[i].trim();
			//System.out.println("line "+i+": "+lines[i]);
			if(!lines[i].startsWith( "#" )){
				String[] temp = lines[i].split( "\t" );
				nameMap.put( temp[0], temp[1] );
				if(temp[0].endsWith( " (Pseudo)" )){
					nameMap.put( temp[0].substring( 0, temp[0].length()-9 ), temp[1] );
				}
			}
		}
		
		
		if(input.getFileContents().getExtension().equalsIgnoreCase( "gff" ) || input.getFileContents().getExtension().equalsIgnoreCase( "gff3" )){
			
			StringBuffer res = new StringBuffer();
			lines = input.getFileContents().getContent().split( "\n" );
			for(int i=0;i<lines.length;i++){
				String[] parts = lines[i].split( "\t" );
				String[] ann = parts[8].split( "\\;" );
				//System.out.println("ann: "+Arrays.toString( ann ));
				for(int j=0;j<ann.length;j++){
					String[] keyval = ann[j].trim().split( "=" );
					//System.out.println("keyval: "+Arrays.toString( keyval ));
					if(keyval[0].equalsIgnoreCase( "id" ) || keyval[0].equalsIgnoreCase( "parent" )){
						if(nameMap.containsKey( keyval[1] )){
							//System.out.println("contains "+keyval[1]);
							keyval[1] = nameMap.get( keyval[1] );
							//System.out.println("set "+keyval[1]);
							ann[j] = keyval[0]+"="+keyval[1];
						}
					}
					
				}
				StringBuffer temp = new StringBuffer();
				for(int j=0;j<ann.length;j++){
					temp.append(ann[j]);
					if(j<ann.length-1){
						temp.append( ";");
					}
				}
				parts[8] = temp.toString();
				for(int j=0;j<parts.length;j++){
					res.append( parts[j] );
					if(j < parts.length-1){
						res.append( "\t" );
					}
				}
				res.append( "\n" );
			}
			
			
			TextResult fres = new TextResult( "Renamed GFF: TALE loci", "Renamed TALEs in GFF format", new FileRepresentation( "", res.toString() ), "gff3", "Rename", null, true );	
			
			protocol.append("Renamed TALEs in GFF file.\n");
			
			return new ToolResult("Result of "+getToolName()+" (GFF)", getToolName()+" on \""+input.getFileContents().getFilename()+"\"", null, new ResultSet(fres), parameters, getToolName(), new Date(System.currentTimeMillis()) );
			
		}else{
			
			ByteArrayOutputStream os = new ByteArrayOutputStream(); 
			
			SimpleNamespace ns = new SimpleNamespace("biojava");
			
			Term gene = RichObjectFactory.getDefaultOntology().getOrCreateTerm( "gene" );
			
			BufferedReader br = new BufferedReader( new StringReader( input.getFileContents().getContent() ) );
			RichSequenceIterator it = RichSequence.IOTools.readGenbankDNA( br, ns );
			LinkedList<RichSequence> li = new LinkedList<RichSequence>();
			while(it.hasNext()){
				RichSequence seq = it.nextRichSequence(); 
				
				Iterator<Feature> featit = seq.features();
				while(featit.hasNext()){
					Feature f = featit.next();
					if(f instanceof RichFeature){
						Iterator<Note> notit = ( (RichFeature)f ).getNoteSet().iterator();
						while(notit.hasNext()){
							Note n = notit.next();
							if(n.getTerm().equals( gene ) && nameMap.containsKey( n.getValue() ) ){
								n.setValue( nameMap.get( n.getValue() ) );
							}
						}
					}
				}
				
				IOTools.writeGenbank( os, seq, ns );
				
			}
			
			String cont = os.toString( "UTF-8" );
			
			TextResult fres2 = new TextResult( "Renamed Genbank: TALEs", "Renamed TALEs in Genbank format", new FileRepresentation( "", cont ), "gb", "Rename", null, true );

			protocol.append("Renamed TALEs in Genbank file.\n");
			
			return new ToolResult("Result of "+getToolName()+" (Genbank)", getToolName()+" on \""+input.getFileContents().getFilename()+"\"", null, new ResultSet(fres2), parameters, getToolName(), new Date(System.currentTimeMillis()) );
			
			
		}
		
		
	}

	@Override
	public String getToolName() {
		return "Rename TALEs in File";
	}

	@Override
	public String getShortName() {
		return "rename";
	}

	@Override
	public String getDescription() {
		return "Renames TALEs in in Genbank or GFF3 file";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( TALEAnalysisTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/RenameTool.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}
	
	@Override
	public String getToolVersion() {
		return "1.1";
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
