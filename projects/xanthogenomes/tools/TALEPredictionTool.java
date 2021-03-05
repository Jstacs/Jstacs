package projects.xanthogenomes.tools;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.SequenceIterator;
import org.biojavax.RichObjectFactory;
import org.biojavax.SimpleNamespace;
import org.biojavax.SimpleNote;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequence.IOTools;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.biojavax.bio.seq.SimplePosition;
import org.biojavax.bio.seq.SimpleRichFeature;
import org.biojavax.bio.seq.SimpleRichLocation;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.bioJava.BioJavaAdapter;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
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
import de.jstacs.utils.Pair;
import projects.xanthogenomes.NHMMer;
import projects.xanthogenomes.Tools;
import projects.xanthogenomes.Tools.Translator;


public class TALEPredictionTool implements JstacsTool {

	public TALEPredictionTool() {
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			FileParameter input = new FileParameter( "Genome", "The input Xanthomonas genome in FastA or Genbank format", "fasta,fa,fas,fna,gb,gbk,genbank", true );

			SimpleParameter strain = new SimpleParameter(DataType.STRING, "Strain", "The name of the strain, will be used for annotated TALEs", false);
			
			SimpleParameter sens = new SimpleParameter(DataType.BOOLEAN, "Sensitive", "Sensitive scan", true,false);
			ToolParameterSet ps  = new ToolParameterSet( getShortName(),  input, strain,sens );
			return ps;
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		
		progress.setLast( 1 );
		progress.setCurrent( 0.0 );
		
		SimpleParameter strainp = (SimpleParameter) parameters.getParameterAt(1);
		String strain = "";
		String strainstr = "";
		if(strainp.isSet()){
			String temp = (String) strainp.getValue();
			if(temp != null && temp.trim().length() > 0){
				strain = " ("+temp+")";
				strainstr = temp+"-";
			}
		}
		
		FileParameter fp = (FileParameter)parameters.getParameterAt( 0 );
		FileRepresentation fr = fp.getFileContents();
		
		boolean sensitive = (boolean) parameters.getParameterAt(2).getValue();

		DataSet ds = null;

		String content = fr.getContent();
		
		RichSequence[] seqs = null;
		SimpleNamespace ns = new SimpleNamespace("biojava");
		
		protocol.append( "Loading input genome...\n" );
		
		AlphabetContainer con = new AlphabetContainer(new DiscreteAlphabet( true, "A", "C", "G", "T", "N", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V" ));
		
		boolean fasta = content.trim().startsWith(">");
		
		if(!fasta){
			try{
				BufferedReader br = new BufferedReader( new StringReader( content ) );

				RichSequenceIterator it = RichSequence.IOTools.readGenbankDNA( br, ns );
				ds = BioJavaAdapter.sequenceIteratorToDataSet( it, null, con );
				br.close();

				br = new BufferedReader( new StringReader( content ) );
				it = RichSequence.IOTools.readGenbankDNA( br, ns );
				LinkedList<RichSequence> li = new LinkedList<RichSequence>();
				while(it.hasNext()){
					RichSequence seq = it.nextRichSequence(); 
					li.add( seq );

				}
				seqs = li.toArray( new RichSequence[0] );
			}catch(Exception e){
				protocol.appendWarning( "... loading failed.\n\n" );
				protocol.appendThrowable(e);
				throw new Exception( "Input not in expected format" );
			}
		}else{
			//e.printStackTrace( );
			try {
				BufferedReader br = new BufferedReader( new StringReader( content ) );
				
				ds = new DataSet( con, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ) );

				SequenceIterator it = BioJavaAdapter.dataSetToSequenceIterator( ds, false, true );
				seqs = new RichSequence[ds.getNumberOfElements()];
				for(int i=0;i<seqs.length;i++){
					seqs[i] = (RichSequence)it.nextSequence();
				}
				
			} catch ( Exception ex ) {
				protocol.appendWarning( "... loading failed.\n\n" );
				protocol.appendThrowable(ex);
				throw new Exception( "Input not in expected format" );
			}	
		}
		
		Pair<DataSet,ArrayList<Integer>[]> pair = preprocess(ds);
		ds = pair.getFirstElement();
		ArrayList<Integer>[] poss = pair.getSecondElement();
		
		protocol.append( "Scanning genome for TALEs...\n" );
		int[][] regions = NHMMer.run( new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/data/repeats.hmm" ) ) , 
				new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream("projects/xanthogenomes/data/starts.hmm")), 
				new InputStreamReader( TALEPredictionTool.class.getClassLoader().getResourceAsStream("projects/xanthogenomes/data/ends.hmm")), ds, progress, sensitive );
		protocol.append( "...finished.\n\n" );
		
		for(int i=0;i<regions.length;i++){
			int seqIdx = regions[i][0];
			int start = regions[i][1];
			int end = regions[i][2];
			boolean contN = false;
			for(int j=0;j<poss[seqIdx].size();j++){
				int idx = poss[seqIdx].get(j);
				//System.out.println(idx+" "+start+" "+end);
				if(idx >= start && idx <= end){
					contN = true;
				}
			}
			if(contN){
				protocol.appendWarning(strainstr+"tempTALE"+(i+1)+" contained \"N\"s in predicted CDS, which have been replaced. Please use with care.\n");
			}
		}
		
		
		protocol.append( "Writing GFF output.\n" );
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<regions.length;i++){
			Sequence seq = ds.getElementAt( regions[i][0] );
			
			String id = null;
			SequenceAnnotation ann = seq.getSequenceAnnotationByType( "unparsed comment line", 0 );

			if(ann != null){
				id = ann.getResultAt( 0 ).getValue().toString().trim();
				/*if(id.indexOf( " " ) > 0){
					id = id.substring( 0, id.indexOf( " " ) );
				}*/
			}else{
				SequenceAnnotation ann2 = seq.getSequenceAnnotationByTypeAndIdentifier( "BioJava RichSequence Annotation", BioJavaAdapter.ANNOTATION_ID );

				if(ann2 != null){
					Result res = ann2.getResultForName( "Name" );
					if(res != null){
						id = res.getValue().toString();
					}
				}
			}
			
			sb.append( id+"\tTALE-prediction\tmRNA\t"+(regions[i][4]+1)+"\t"+regions[i][5]+"\t.\t"+(regions[i][3] < 0 ? "-" : "+")+"\t.\tId="+strainstr+"tempTALE"+(i+1)+(regions[i][6] == 0 ? "" : "; Note=putative pseudo gene")+"\n" );
			sb.append( id+"\tTALE-prediction\tCDS\t"+(regions[i][1]+1)+"\t"+regions[i][2]+"\t.\t"+(regions[i][3] < 0 ? "-" : "+")+"\t.\tParent="+strainstr+"tempTALE"+(i+1)+"\n" );
			
		}
				
		TextResult fres = new TextResult( "GFF: TALE predictions"+strain, "TALE predictions in GFF format", new FileRepresentation( "", sb.toString() ), "gff3", "TALE Prediction", null, true );
		
		
		protocol.append( "Writing Genbank output.\n" );
		for(int i=0;i<regions.length;i++){
			RichFeature.Template temp = new RichFeature.Template();
			temp.location = new SimpleRichLocation( new SimplePosition( regions[i][1]+1), new SimplePosition( regions[i][2] ), i, regions[i][3] < 0 ? RichLocation.Strand.NEGATIVE_STRAND : RichLocation.Strand.POSITIVE_STRAND );
			temp.source = "TALE-prediction";
			temp.type = "CDS";
			temp.annotation = new SimpleAnnotation();
			temp.featureRelationshipSet = new HashSet();
			temp.rankedCrossRefs = new HashSet();
			SimpleRichFeature feat = new SimpleRichFeature( seqs[regions[i][0]], temp );
			
			feat.getNoteSet().add( new SimpleNote( RichObjectFactory.getDefaultOntology().getOrCreateTerm( "gene" ), strainstr+"tempTALE"+(i+1), 1 ) );
		
			seqs[regions[i][0]].getFeatureSet().add( feat );
			
			temp.type = "mRNA";
			temp.location = new SimpleRichLocation( new SimplePosition( regions[i][4]+1), new SimplePosition( regions[i][5] ), i, regions[i][3] < 0 ? RichLocation.Strand.NEGATIVE_STRAND : RichLocation.Strand.POSITIVE_STRAND );
			feat = new SimpleRichFeature( seqs[regions[i][0]], temp );
			feat.getNoteSet().add( new SimpleNote( RichObjectFactory.getDefaultOntology().getOrCreateTerm( "gene" ), strainstr+"tempTALE"+(i+1), 1 ) );
			if(regions[i][6] == 1){
				feat.getNoteSet().add( new SimpleNote( RichObjectFactory.getDefaultOntology().getOrCreateTerm( "note" ), "putative pseudo gene", 2 ) );
			}
			
			seqs[regions[i][0]].getFeatureSet().add( feat );
		}
		
		ByteArrayOutputStream os = new ByteArrayOutputStream(); 
		
		for(int i=0;i<seqs.length;i++){
			IOTools.writeGenbank( os, seqs[i], ns );
			
		}
		
		String cont = os.toString( "UTF-8" );
		
		TextResult fres2 = new TextResult( "Genbank: TALE predictions"+strain, "Annotated TALEs in Genbank format", new FileRepresentation( "", cont ), "gb", "TALE Prediction", null, true );
		
		StringBuffer pseudo = new StringBuffer();
		StringBuffer dna = new StringBuffer();
		StringBuffer prot = new StringBuffer();
		protocol.append( "Writing FastA outputs.\n" );
		for(int i=0;i<regions.length;i++){
			//System.out.println(i+" "+regions[i][1]+" "+regions[i][2]);
			Sequence seq = ds.getElementAt( regions[i][0] ).getSubSequence( regions[i][1], regions[i][2]-regions[i][1] );
			if(regions[i][3] < 0){
				seq = seq.reverseComplement();
			}
			String posString = "["+regions[i][1]+"-"+regions[i][2]+":"+regions[i][3]+"]";
			dna.append( ">"+strainstr+"tempTALE"+(i+1)+(regions[i][6] == 1 ? " (Pseudo)" : "")+" "+posString+"\n"+seq+"\n" );
			Sequence seq2 = Translator.DEFAULT.translate( seq, 0 );
			prot.append( ">"+strainstr+"tempTALE"+(i+1)+(regions[i][6] == 1 ? " (Pseudo)" : "")+" "+posString+"\n"+seq2+"\n" );
			
			if(regions[i][6] == 1){
				Sequence pseudoDNA = ds.getElementAt( regions[i][0] ).getSubSequence( regions[i][4], regions[i][5]-regions[i][4] );
				if(regions[i][3] < 0){
					pseudoDNA = pseudoDNA.reverseComplement();
				}
				posString = "["+regions[i][4]+"-"+regions[i][5]+":"+regions[i][3]+"]";
				pseudo.append( ">"+strainstr+"tempTALE"+(i+1)+(regions[i][6] == 1 ? " (Pseudo)" : "")+" "+posString+"\n"+pseudoDNA+"\n" );
				for(int j=0;j<3;j++){
					Sequence pseudoProt = Tools.Translator.DEFAULT.translate(pseudoDNA, j);
					pseudo.append( ">"+strainstr+"tempTALE"+(i+1)+(regions[i][6] == 1 ? " (Pseudo)" : "")+" frame: "+j+" "+posString+"\n"+pseudoProt+"\n" );
				}
			}
		}
		
		TextResult fres3 = new TextResult( "TALE DNA sequences"+strain, "The DNA sequences of the TALE CDS", new FileRepresentation( "", dna.toString() ), "fasta", "TALE Prediction", "fasta/dna", true );
		TextResult fres4 = new TextResult( "TALE protein sequences"+strain, "The protein sequences of the TALE CDS", new FileRepresentation( "", prot.toString() ), "fasta", "TALE Prediction", "fasta/as", true );
		
		ResultSet set = null;
		if(pseudo.length() > 0){
			TextResult fres5 = new TextResult( "TALE pseudo gene matches"+strain, "The complete matching sequences of the TALE pseudo genes as DNA and translated in all three reading frames", new FileRepresentation( "", pseudo.toString() ), "fasta", "TALE Prediction", "fasta/as", true );
			set = new ResultSet( new Result[]{fres,fres2,fres3, fres4, fres5} );
		}else{
			set = new ResultSet( new Result[]{fres,fres2,fres3, fres4} );
		}
		
		
		String file = "";
		if(fr.getFilename() != null){
			file = " on "+fr.getFilename();
		}
		
		return new ToolResult("Result of "+getToolName()+strain, getToolName()+file, null, set, parameters, getToolName(), new Date(System.currentTimeMillis()));
		
	}

	private Pair<DataSet,ArrayList<Integer>[]> preprocess(DataSet ds) throws IllegalArgumentException, WrongAlphabetException, EmptyDataSetException {
		Sequence[] seqs = ds.getAllElements();
		Pattern pat = Pattern.compile("[^ACGTacgt]");
		ArrayList<Integer>[] poss = new ArrayList[ds.getNumberOfElements()];
		for(int i=0;i<seqs.length;i++){
			poss[i] = new ArrayList<Integer>();
			String seqstr = seqs[i].toString();
			SequenceAnnotation[] anns = seqs[i].getAnnotation();
			Matcher m = pat.matcher(seqstr);
			while(m.find()){
				int pos = m.start();
				poss[i].add(pos);
			}
			seqstr = m.replaceAll("A");
			seqs[i] = Sequence.create(DNAAlphabetContainer.SINGLETON, anns, seqstr, "");
		}
		return new Pair<DataSet, ArrayList<Integer>[]>( new DataSet(ds.getAnnotation(),seqs), poss);
	}

	@Override
	public String getToolName() {
		return "TALE Prediction";
	}

	@Override
	public String getShortName() {
		return "predict";
	}

	@Override
	public String getDescription() {
		return "Predicts TALE in a genome";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( TALEPredictionTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/TALEPredictionTool.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}
	
	@Override
	public String getToolVersion() {
		return "1.4.2";
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
