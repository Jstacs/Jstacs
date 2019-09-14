package projects.xanthogenomes.tools;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

import de.jstacs.DataType;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
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
import projects.tals.ScanForTBSCLI;
import projects.tals.TALgetterDiffSM;
import projects.xanthogenomes.FamilyGroupPlotter;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALEFamilyBuilder;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;


public class LoadAndViewClassesTool implements JstacsTool {

	public LoadAndViewClassesTool() {
	}

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter builderFile = new FileParameter( "Class builder", "TALE class builder definition", "xml", true );
		builderFile.setExtendedType( TALEFamilyBuilder.class.getName() );
		
		try {
			SelectionParameter sp = new SelectionParameter( DataType.PARAMETERSET, new String[]{"Download current definition","Load from local file"}, new Object[]{new SimpleParameterSet( ), new SimpleParameterSet( builderFile )}, "Class definition source", "Download current class definition (requires internet connection) or load definition from local file", true );
			
			return new ToolParameterSet( getShortName(),  sp );
		} catch ( Exception e ) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		
	}

	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		
		progress.setLast( 2.0 );
		progress.setCurrent( 0.0 );
		
		SelectionParameter sp = (SelectionParameter)parameters.getParameterAt( 0 );
		
		TALEFamilyBuilder builder = null;
		
		
		TALgetterDiffSM model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLI.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/talfinder_obg2_hyp_bg.xml" ) ), "model" );
		
		
		
		LinkedList<Result> ress = new LinkedList<Result>();
		
		if(sp.getSelected() == 0){
			protocol.append( "Downloading current class definition...\n" );
		//	Scanner scan = new Scanner(new URL("http://www.jstacs.de/downloads/class_definitions_current.xml.gz").openStream(),"UTF-8");//TODO URL
		//	String text = scan.useDelimiter( "\\A" ).next();
		//	scan.close();
			InputStream is = new URL("http://www.jstacs.de/downloads/class_definitions_current2.xml.gz").openStream();
			GZIPInputStream gzis = new GZIPInputStream(is);
			Scanner scan = new Scanner(gzis);
			String text = scan.useDelimiter( "\\A" ).next();
			scan.close();
			
			protocol.append( "...finished.\n\nRe-building class builder...\n" );
			builder = new TALEFamilyBuilder( new StringBuffer(text) );
			ress.add( new TextResult( "Class builder download", "TALE class builder definition", new FileRepresentation( "", text ), "xml", "Downloaded", TALEFamilyBuilder.class.getName(), true ) );
			protocol.append("...finished.\n\n");
		}else{
			
			ParameterSet ps = (ParameterSet)sp.getValue();
			protocol.append( "Loading class builder...\n" );
			builder = new TALEFamilyBuilder(new StringBuffer( ((FileParameter)ps.getParameterAt( 0 )).getFileContents().getContent()));
			protocol.append("...finished.\n\n");
		}
		progress.setCurrent( 1.0 );
		
		builder.setToOld();
			
		
		ClusterTree<TALEFamily> famTree = builder.clusterFamilies();
		TALEFamily[] fams = builder.getFamilies();
		Arrays.sort( fams );
			
		
		TALE[] tales = builder.getAllTALEs();
		LinkedList<Result> reports = new LinkedList<Result>();
		StringBuffer allTales = new StringBuffer();
		for(int i=0;i<tales.length;i++){
			allTales.append(">"+tales[i].getId()+" "+tales[i].annotationToString()+"\n");
			allTales.append(tales[i].getStart());
			for(int j=0;j<tales[i].getNumberOfRepeats();j++){
				allTales.append(tales[i].getRepeat(j).getRepeat());
			}
			allTales.append(tales[i].getEnd());
			allTales.append("\n");
		}
		reports.add(new TextResult("All TALE protein sequences", "Sequences of all TALEs in the class builder", new FileRepresentation("", allTales.toString()), "fasta", "LoadAndViewClasses", "fasta/as", true));
		
		HashMap<String,StringBuffer> talesByStrain = new HashMap<String, StringBuffer>();
		allTales = new StringBuffer();
		for(int i=0;i<tales.length;i++){
			TALE dna = tales[i].getDnaOriginal();
			if(dna != null){
				StringBuffer temp = new StringBuffer();
				temp.append(">"+dna.getId()+" "+dna.annotationToString()+"\n");
				temp.append(dna.getStart());
				for(int j=0;j<dna.getNumberOfRepeats();j++){
					temp.append(dna.getRepeat(j).getRepeat());
				}
				temp.append(dna.getEnd());
				temp.append("\n");
				allTales.append(temp);
				
				String strain = dna.getStrain()+"";
				if(!talesByStrain.containsKey(strain)){
					talesByStrain.put(strain, new StringBuffer());
				}
				talesByStrain.get(strain).append(temp);
			}
			
		}
		reports.add(new TextResult("All TALE DNA sequences", "Sequences of all TALEs in the class builder", new FileRepresentation("", allTales.toString()), "fasta", "LoadAndViewClasses", "fasta/dna", true));

//System.out.println("by strain");
		LinkedList<Result> byStrainList = new LinkedList<Result>();
		Iterator<String> it = talesByStrain.keySet().iterator();
		while(it.hasNext()){
			String strain = it.next();
			String strainn = "null".equals(strain) ? "unknown" : strain;
			byStrainList.add(new TextResult("TALE DNA sequences ("+strainn+")", "Sequences of all TALEs for strain "+strainn, new FileRepresentation("", talesByStrain.get(strain).toString()), "fasta", "LoadAndViewClasses", "fasta/dna", true));
		}
		
		reports.add(new ResultSetResult("TALEs by strain", "The TALE sequences of all strains", null, new ResultSet(byStrainList)));
		
		
		/*StringBuffer taleList = new StringBuffer();
		for(int i=0;i<tales.length;i++){
			taleList.append(tales[i].getId()+"\t"+tales[i].annotationToColumns()+"\n");
		}
		reports.add(new TextResult("List of TALEs", "The list of all TALEs in the class builder", new FileRepresentation("", taleList.toString()), "txt", "LoadAndViewClasses", "txt", true));
		*/
//System.out.println("List 1");	
		LinkedList<ResultSet> taleList = new LinkedList<ResultSet>();
		HashMap<String,String[]> strainMap = new HashMap<String, String[]>();
		LinkedList<ResultSet> strainRes = new LinkedList<ResultSet>();
		int woStrain = 0;
		for(int i=0;i<tales.length;i++){
			taleList.add(tales[i].annotationToResultSet());
			
			String strain = tales[i].getStrain();
			if(strain != null){
				String acc = tales[i].getAccession();
				if(acc == null){
					acc = "";
				}
				if(strainMap.containsKey(acc)){
					String[] cont = strainMap.get(acc);
					if(!cont[1].equals(acc)){
						strainRes.add(new ResultSet(new Result[]{new CategoricalResult("Strain", "", strain),
								new CategoricalResult("Accession", "", acc)}));
					}
				}else{
					strainRes.add(new ResultSet(new Result[]{new CategoricalResult("Strain", "", strain),
							new CategoricalResult("Accession", "", acc)}));
					strainMap.put(acc, new String[]{strain,acc});
				}
			}else{
				woStrain++;
			}
		}
		if(woStrain > 0){
			strainRes.add(new ResultSet(new Result[]{new CategoricalResult("Strain", "", "No strain annotation"),
					new CategoricalResult("Accession", "", woStrain+" TALEs")}));
		}
		ListResult lr = new ListResult("List of TALEs", "The list of all TALEs in the class builder", null, taleList.toArray(new ResultSet[0]));
		reports.add(lr);
		
		ListResult lr2 = new ListResult("List of strains", "The list of all strains in the class builder", null, strainRes.toArray(new ResultSet[0]));
		reports.add(lr2);

		StringBuffer familyReport = new StringBuffer();
		for(int i=0;i<fams.length;i++){
			/*familyReport.append("Class "+fams[i].getFamilyId()+":\n\n");
			familyReport.append("Class tree:\n"+fams[i].getTree().toNewick()+"\n\n");
			familyReport.append("Class members:\n\n");
			TALE[] mems = fams[i].getFamilyMembers();
			for(int j=0;j<mems.length;j++){
				familyReport.append(mems[j].getId()+"\t"+mems[j].getRvdSequence()+"\n");
			}
			familyReport.append( fams[i].inducedMultipleAlignmentToString()+"\n");*/
			familyReport.append(fams[i].toString(model, builder));
			familyReport.append("\n\n#####################################################\n\n");
		}

		reports.add(new TextResult("List of classes", "The list of all TALEs classes and corresponding distance trees", new FileRepresentation("", familyReport.toString()), "txt", "LoadAndViewClasses", "txt", true));

		ResultSet repSet = new ResultSet(reports.toArray(new Result[0]));

		ress.add(new ResultSetResult("Lists of classes, strains and TALEs", "Details on the classes, strains and TALEs in the current class builder", null, repSet));
		
		

		ress.add( new PlotGeneratorResult( "Tree of classes", "The tree of class similarities", new FamilyGroupPlotter.FamilyGroupPlotGenerator( famTree ), true ) ); 

		
		LinkedList<Result> classTrees = new LinkedList<Result>();
		
		protocol.append( "Creating class plots for...\n" );
		for(int i=0;i<fams.length;i++){
			protocol.append( fams[i].getFamilyId()+"\n" );
			PlotGeneratorResult pgr = new PlotGeneratorResult( "Class "+fams[i].getFamilyId(), "Plot of the tree of the TALEs in this class", fams[i], true );
			classTrees.add( pgr );
			progress.setCurrent( 1.0 + i/(double)fams.length );
		}
		
		ress.add(new ResultSetResult("Class trees", "Trees of all TALE classes", null, new ResultSet(classTrees)));
		
		ResultSet set = new ResultSet( ress.toArray( new Result[0] ) );
		
		return new ToolResult("Result of "+getToolName(), "", null, set, parameters, getToolName(), new Date(System.currentTimeMillis())); 
		
	}

	@Override
	public String getToolName() {
		return "Load and View TALE Classes";
	}

	@Override
	public String getShortName() {
		return "loadAndView";
	}

	@Override
	public String getDescription() {
		return "Downloads and displays current TALE class definition";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( LoadAndViewClassesTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/LoadAndViewClassesTool.txt" ) ).toString();
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
		return null;
	}


}
