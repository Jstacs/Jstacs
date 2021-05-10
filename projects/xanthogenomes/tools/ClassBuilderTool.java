package projects.xanthogenomes.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.data.DataSet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.PlotGeneratorResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.Pair;
import projects.tals.ScanForTBSCLI;
import projects.tals.TALgetterDiffSM;
import projects.xanthogenomes.BuildFamilies;
import projects.xanthogenomes.FamilyGroupPlotter;
import projects.xanthogenomes.RVDAlphabetContainer;
import projects.xanthogenomes.SplitTALEs;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALEFamilyBuilder;
import projects.xanthogenomes.BuildFamilies.FamilyResult;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;
import projects.xanthogenomes.Tools.ProteinAlphabetContainer;


public class ClassBuilderTool implements JstacsTool {

	public ClassBuilderTool() {
		
	}

	@Override
	public ToolParameterSet getToolParameters() {
		
		try{
			FileParameter fp = new FileParameter("TALE sequences","The sequences of the TALEs (DNA or protein), or \"TALE DNA parts\" or \"TALE Protein parts\" output of \"TALE Analysis\", or RVD sequences.","fasta,fa,fas",true);
			fp.setExtendedType( "fasta/dna" );

			SimpleParameter cut = new SimpleParameter( DataType.DOUBLE, "Cutoff", "Cutoff value defining the maximum distance of a TALE class", true, new NumberValidator<Comparable<Double>>( 0.0, Double.MAX_VALUE ), 5.0 );

			SimpleParameter pval = new SimpleParameter( DataType.DOUBLE, "Significance level", "Cutoff value on the p-value representing alignment significance", true, new NumberValidator<Comparable<Double>>( 0.0, 1.0 ), 0.01 );

			return new ToolParameterSet( getShortName(), fp,cut,pval);
			
		}catch(ParameterException e){
			e.printStackTrace();
			return null;
		}
	}

	public static TALE[] readProteinTALEs(FileRepresentation fr, Protocol protocol) throws Exception{
		
		String content = fr.getContent();
		
		try{
			TALE[] res = TALEAnalysisTool.parseTALEsFromParts( content, protocol );
			protocol.append( "Loaded TALEs from \"TALE Analysis\" parts.\n" );
			return res;
		}catch(Exception e){
			
			//e.printStackTrace( );
			
			try{
				DataSet ds = null;
				try{
					BufferedReader br = new BufferedReader( new StringReader( content ) );
					ds = new DataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ) );
				}catch(Exception ex){
					BufferedReader br = new BufferedReader( new StringReader( content ) );
					ds = new DataSet( ProteinAlphabetContainer.SINGLETON, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ) );
				}

				LinkedList<TALE> talelist = new LinkedList<TALE>();

				HashSet<String> ids = new HashSet<String>();

				for(int i=0;i<ds.getNumberOfElements();i++){
					String id = ds.getElementAt( i ).getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultAt( 0 ).getValue().toString();
					protocol.append( id+"\n" );
					TALE[] tales = SplitTALEs.split( id, ds.getElementAt( i ), protocol );

					TALE prot = tales[1];
					if(prot != null){
						int c = 0;
						String tempId = prot.getId();
						while(ids.contains( prot.getId() )){
							if(c == 0){
								protocol.appendWarning( "Duplicate ID "+tempId+"\n" );
							}
							c++;
							prot.setId( tempId+c );
						}
					}
					if(prot != null && prot.getNumberOfRepeats() > 0){
						talelist.add( prot );
					}else if(prot != null){
						protocol.appendWarning( "Removed putative pseudo gene "+prot.getId()+", because it has zero repeats.\n" );
						//System.err.println( "Removed putative pseudo gene "+prot.getId()+"." );
					}else{
						protocol.appendWarning( "Removed putative pseudo gene "+id+", because it could not be translated.\n" );
					}
				}
				protocol.append( "Loaded TALEs from complete sequences.\n" );

				TALE[] ttales = talelist.toArray( new TALE[0] );

				return ttales;
			}catch(Exception ex){
				BufferedReader br = new BufferedReader( new StringReader( content ) );
				DataSet ds = new DataSet( RVDAlphabetContainer.SINGLETON, new SparseStringExtractor( br, '>', "", new SimpleSequenceAnnotationParser() ), "-" );
				TALE[] tales = new TALE[ds.getNumberOfElements()];
				for(int i=0;i<ds.getNumberOfElements();i++){
					String id = ds.getElementAt( i ).getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultAt( 0 ).getValue().toString();
					protocol.append( id+"\n" );
					tales[i] = new TALE(id, ds.getElementAt(i), true, true);
				}
				protocol.appendWarning("Loaded TALEs from RVD sequences. Some properties of TALEs (aberrant repeats, codon mismatches) cannot be displayed due to lacking information.\n");
				return tales;
			}
		}
	}
	
	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		
		progress.setLast( 1.0 );
		progress.setCurrent( 0.0 );
		
		FileRepresentation fr = ((FileParameter)parameters.getParameterAt( 0 )).getFileContents();
		double cut = (Double)parameters.getParameterAt(1).getValue();
		double pval = (Double)parameters.getParameterAt(2).getValue();
		
		TALE[] ttales = readProteinTALEs( fr, protocol );
		
		progress.setLast( 0.3 );
		
		protocol.append( "Building classes.\n" );
		
		Pair<TALEFamilyBuilder,FamilyResult[]> res = BuildFamilies.build( ttales, cut, pval );
		
		progress.setCurrent( 0.7 );
		
		FamilyResult[] famRes = res.getSecondElement();
		Arrays.sort( famRes );
		
		TALgetterDiffSM model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLI.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/talfinder_obg2_hyp_bg.xml" ) ), "model" );
		
		ClusterTree<TALEFamily> famTree = res.getFirstElement().clusterFamilies();
		
		progress.setCurrent( 0.9 );
		
		Result[] ress = new Result[famRes.length+2];
		
		ress[0] = new TextResult( "Class builder", "TALE class builder definition", new FileRepresentation( "", res.getFirstElement().toXML().toString() ), "xml", "TALE Class Builder", TALEFamilyBuilder.class.getName(), true );
		ress[1] = new PlotGeneratorResult( "Tree of classes", "The tree of class similarities", new FamilyGroupPlotter.FamilyGroupPlotGenerator( famTree ), true ); 
		
		protocol.append( "Generating reports and plots for "+famRes.length+" classes...\n" );
		
		for(int i=0;i<famRes.length;i++){
			//famRes[i].setSimilarityGroup( groups.getFirstElement() );
			TALEFamily fam = famRes[i].getFamily();
			protocol.append( fam.getFamilyId()+"\n" );
			
			PlotGeneratorResult pgr = new PlotGeneratorResult( "Class tree for "+fam.getFamilyId(), "Plot of the tree of the TALEs in this class", fam, true );
			
			//fam.plotFamilyToFile(outpath+"/family_"+(i+1), adaptor );
			String report = famRes[i].toString(model, res.getFirstElement());
			
			TextResult fileres = new TextResult( "Class report for "+fam.getFamilyId(), "Report for class "+fam.getFamilyId(), new FileRepresentation( "", report ), "txt", "TALE Class Builder", null, false );
			
			ResultSetResult rsr = new ResultSetResult( "Class "+fam.getFamilyId(), "Collection of results for class "+fam.getFamilyId(), null, new ResultSet( new Result[]{fileres,pgr} ) );
			
			ress[i+2] = rsr;
			
			//double[][] specs = fam.getSpecificityProfile( model );
			//int w = SeqLogoPlotter.getWidth( 300, specs );
			/*Graphics2D graph = adaptor.getGraphics( w, 300 );
			String[] labels = new String[specs.length];
			for(int j=0;j<labels.length;j++){
				labels[j] = j+"";
			}
			SeqLogoPlotter.plotLogo( graph, 300, specs, labels, "Position", "bits" );
			adaptor.generateOutput( outpath+"/family_"+(i+1)+"_specificity."+adaptor.getGraphicsExtension() );*/
		}
		
		ResultSet set = new ResultSet( ress );
		return new ToolResult("Result of "+getToolName(), "", null, set, parameters, getToolName(), new Date(System.currentTimeMillis()));
		
	}

	@Override
	public String getToolName() {
		return "TALE Class Builder";
	}

	@Override
	public String getShortName() {
		return "build";
	}

	@Override
	public String getDescription() {
		return "Creates classes from a set of input TALEs";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( ClassBuilderTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/ClassBuilderTool.txt" ) ).toString();
		} catch ( IOException e ) {
			e.printStackTrace();
			return "";
		}
	}

	@Override
	public String getToolVersion() {
		return "1.4.1";
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
		return new String[] {"@article{grau16annotale,\n" + 
				"	title = {{AnnoTALE}: bioinformatics tools for identification, annotation, and nomenclature of {TALEs} from \\emph{Xanthomonas} genomic sequences},\n" + 
				"	author = {Grau, Jan and Reschke, Maik and Erkes, Annett and Streubel, Jana and Morgan, Richard D. and Wilson, Geoffrey G. and Koebnik, Ralf and Boch, Jens},\n" + 
				"	journal = {Scientific Reports},\n" + 
				"	year = {2016},\n" + 
				"	volume = {6},\n" + 
				"	pages = {21077},\n" + 
				"	doi = {10.1038/srep21077}\n" + 
				"	}\n"};
	}

	
}
