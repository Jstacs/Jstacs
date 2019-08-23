package projects.xanthogenomes.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Date;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
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
import projects.xanthogenomes.SplitTALEs;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALE.Repeat;
import projects.xanthogenomes.TALE.Type;
import projects.xanthogenomes.Tools.ProteinAlphabetContainer;
import projects.xanthogenomes.Tools.Translator;


public class TALEAnalysisTool implements JstacsTool {

	public TALEAnalysisTool() {
		
	}

	@Override
	public ToolParameterSet getToolParameters() {

		try{
			SimpleParameter name = new SimpleParameter(DataType.STRING, "Name", "A name for this run of "+getToolName(), false);
			FileParameter fp = new FileParameter("TALE DNA sequences","The DNA sequences of the TALEs","fasta,fa,fas",true);
			fp.setExtendedType( "fasta/dna" );

			return new ToolParameterSet( getShortName(),  name, fp );
		}catch(Exception e){
			return null;
		}
	}

	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		String name = (String) parameters.getParameterAt(0).getValue();
		FileParameter fp = (FileParameter)parameters.getParameterAt( 1 );
		FileRepresentation fr = fp.getFileContents();
		
		BufferedReader br = new BufferedReader( new StringReader( fr.getContent() ) );
		DataSet ds = new DataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ) );
		
		
		StringBuffer dnaTALEs = new StringBuffer();
		StringBuffer protTALEs = new StringBuffer();
		StringBuffer rvdTALEs = new StringBuffer();
		
		progress.setLast( ds.getNumberOfElements() );
		progress.setCurrent( 0 );
		
		protocol.append( "Analyzing TALEs:\n" );
		for(int i=0;i<ds.getNumberOfElements();i++){
			String id = ds.getElementAt( i ).getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultAt( 0 ).getValue().toString();
			protocol.append( id+"\n" );
			TALE[] tales = SplitTALEs.split( id, ds.getElementAt( i ), protocol );
			TALE dna = tales[0];
			
			
			if(dna != null){
				dnaTALEs.append( ">"+id+": N-terminus\n"+dna.getStart().toString()+"\n" );
				for(int j=0;j<dna.getNumberOfRepeats();j++){
					dnaTALEs.append( ">"+id+": repeat "+(j+1)+"\n"+dna.getRepeat( j ).getRepeat().toString()+"\n" );
				}
				dnaTALEs.append( ">"+id+": C-terminus\n"+dna.getEnd().toString()+"\n" );
			}else{
				protocol.appendWarning("TALE "+id+" could not be analyzed, splitting into regions failed.\n");
			}
			
			TALE prot = tales[1];
			if(prot != null){
				Sequence rvds = tales[1].getRvdSequence();

				protTALEs.append( ">"+id+": N-terminus\n"+prot.getStart().toString()+"\n" );
				for(int j=0;j<prot.getNumberOfRepeats();j++){
					protTALEs.append( ">"+id+": repeat "+(j+1)+"\n"+prot.getRepeat( j ).getRepeat().toString()+"\n" );
				}
				protTALEs.append( ">"+id+": C-terminus\n"+prot.getEnd().toString()+"\n" );

				String rvdStr = rvds.toString( "-", 0, rvds.getLength() );
				
				if(prot.containsAberrantRepeat()){
					String[] rvdAr = rvdStr.split("-");
					for(int j=0;j<prot.getNumberOfRepeats()-1;j++){
						if(prot.getRepeat(j).getType() == Type.LONG || prot.getRepeat(j).getType() == Type.SHORT){
							rvdAr[j] = rvdAr[j].toLowerCase();
						}
					}
					rvdStr = String.join("-", rvdAr);
				}
				
				rvdTALEs.append( ">"+id+"\n"+rvdStr+"\n" );
				progress.setCurrent( i );
			}else{
				protocol.appendWarning("Protein version of TALE "+id+" could not be analyzed, splitting into regions failed.\n");
			}
		}
		
		protocol.append( "\nWriting outputs.\n" );
		
		String names = (name == null || name.length()==0 ? "" : " ("+name+")");
		TextResult fr1 = new TextResult( "TALE DNA parts"+names, "Parts of the TALE (N-/C-terminus, repeats) as DNA sequences", new FileRepresentation( "", dnaTALEs.toString() ), "fasta", "TALE Analysis", "fasta/dna",true );
		TextResult fr2 = new TextResult( "TALE Protein parts"+names, "Parts of the TALE (N-/C-terminus, repeats) as protein sequences", new FileRepresentation( "", protTALEs.toString() ), "fasta", "TALE Analysis", "fasta/as",true );
		TextResult fr3 = new TextResult( "TALE RVDs"+names, "Sequence of RVDs of the TALEs", new FileRepresentation( "", rvdTALEs.toString() ), "fasta", "TALE Analysis", "fasta/rvd",true );
		
		
		ResultSet set = new ResultSet( new Result[]{fr1,fr2,fr3} );
		return new ToolResult("Result of "+getToolName()+names, getToolName()+" on \""+fr.getFilename()+"\"", null, set, parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}
	
	
	public static TALE[] parseTALEsFromParts(String content, Protocol protocol) throws Exception{
		
		DataSet ds = null;
		try{
			BufferedReader br = new BufferedReader( new StringReader( content ) );
			ds = new DataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ) );
		}catch(Exception e){
			//e.printStackTrace( );
			try{
				BufferedReader br = new BufferedReader( new StringReader( content ) );
				ds = new DataSet( ProteinAlphabetContainer.SINGLETON, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ) );
			}catch(Exception ex){
				throw e;
			}
		}
		
		
		String lastid = null;
		Sequence nterm = null;
		Sequence cterm = null;
		LinkedList<Sequence> repeats = new LinkedList<Sequence>();
		LinkedList<TALE> tales = new LinkedList<TALE>();
		for(int i=0;i<ds.getNumberOfElements();i++){
			Sequence seq = ds.getElementAt( i );
			String name = seq.getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultAt( 0 ).getValue().toString();
			int colon = name.lastIndexOf( ":" );
			String id = name.substring( 0, colon ).trim();
			
			if(lastid != null && !lastid.equals( id ) ){
				if(repeats.size() > 0){
					Repeat[] repeatss =new Repeat[repeats.size()];
					for(int j=0;j<repeats.size();j++){
						repeatss[j] = new Repeat( repeats.get( j ) );
					}
					TALE t = new TALE( lastid, nterm, repeatss, cterm, false );
					if(ds.getAlphabetContainer().checkConsistency( DNAAlphabetContainer.SINGLETON )){
						try{
							t = t.getTranslatedTALE( Translator.DEFAULT );
						}catch(Exception e){
							protocol.appendWarning("Could not translate TALE "+t.getId()+". Reason:\n"+e.getMessage()+"\n");
							t = null;
						}
					}
					if(t != null){
						tales.add( t );
					}
				}else{
					protocol.appendWarning("Removed TALE "+lastid+", because it has zero repeats.\n");
				}
				nterm = null;
				cterm = null;
				repeats.clear();
			}
			
			String part = name.substring( colon+1 ).trim();
			if("N-terminus".equals( part )){
				nterm = seq;
			}else if("C-terminus".equals( part )){
				cterm = seq;
			}else if(part.startsWith( "repeat" )){
				repeats.add( seq );
			}else{
				throw new Exception();
			}
			lastid = id;
		}	
		
		if(cterm != null && nterm != null && repeats.size() > 0 ){
			Repeat[] repeatss =new Repeat[repeats.size()];
			for(int j=0;j<repeats.size();j++){
				repeatss[j] = new Repeat( repeats.get( j ) );
			}
			TALE t = new TALE( lastid, nterm, repeatss, cterm, false );
			if(ds.getAlphabetContainer().checkConsistency( DNAAlphabetContainer.SINGLETON )){
				t = t.getTranslatedTALE( Translator.DEFAULT );
			}
			tales.add( t );
		}
		
		
		return tales.toArray( new TALE[0] );
		
	}
	

	@Override
	public String getToolName() {
		return "TALE Analysis";
	}

	@Override
	public String getShortName() {
		return "analyze";
	}

	@Override
	public String getDescription() {
		return "Analyzes TALE structure and RVDs";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( TALEAnalysisTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/TALEAnalysisTool.txt" ) ).toString();
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
