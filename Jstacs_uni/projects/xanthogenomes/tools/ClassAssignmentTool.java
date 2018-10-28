package projects.xanthogenomes.tools;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
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
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.Pair;
import projects.tals.ScanForTBSCLI;
import projects.tals.TALgetterDiffSM;
import projects.xanthogenomes.BuildFamilies;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALEFamilyBuilder;
import projects.xanthogenomes.BuildFamilies.FamilyResult;
import projects.xanthogenomes.TALEFamilyBuilder.FamilyIdGenerator;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;


public class ClassAssignmentTool implements JstacsTool {

	private static LinkedList<String> alreadyGiven = new LinkedList<String>();
	
	private class SchemaFamilyIdGenerator implements FamilyIdGenerator{

		private String orgSuff;
		private String accession;
		
		private String[][] lastNameMap;
		
		
		
		public SchemaFamilyIdGenerator(String orgSuff, String accession){
			this.orgSuff = orgSuff;
			this.accession = accession;
			this.lastNameMap = new String[0][2];
		}
		
		private String[] getFamIDs(){
			String[] temp = new String[26*26];
			for(int i=65,k=0;i<91;i++){
				for(int j=65;j<91;j++,k++){
					temp[k] = "Tal"+(((char)i)+"")+(((char)j)+"");
				}
			}
			return temp;
		}
		
		
		@Override
		public void setFamilyIDs( TALEFamily[] families, TALEFamilyBuilder builder ) {
			
			LinkedList<String[]> nameMap = new LinkedList<String[]>();
			
			HashSet<String> reserved = new HashSet<String>();
			reserved.addAll(alreadyGiven);
			String[] re = builder.getReservedNames();
			if(re != null){
				for(int i=0;i<re.length;i++){
					reserved.add( re[i] );
				}
			}
			
			String regex = "^Tal[A-Z]{2}[0-9]+";
			Pattern pat = Pattern.compile( regex );
			
			for(int i=0;i<families.length;i++){

				String famId = families[i].getFamilyId();
				if(famId != null){
					reserved.add( famId );
					TALE[] mems = families[i].getFamilyMembers();
					for(int j=0;j<mems.length;j++){
						if(!mems[j].isNew()){
							String id = mems[j].getId();
							if(id.matches( regex+".*" )){
								Matcher matcher = pat.matcher( id );
								if(matcher.find()){
									reserved.add( matcher.group() );
								}
							}
						}
					}
				}
			}
			
			String[] syms = getFamIDs();
			int off = 0;
			while( reserved.contains( syms[off] ) ){
				off++;
			}
			
			
			for(int i=0;i<families.length;i++){
				String famId = families[i].getFamilyId();
				if(famId == null){
					while(reserved.contains( syms[off] )){
						off++;
					}
					famId = syms[off];
					families[i].setFamilyId( famId );
					reserved.add( famId );
					alreadyGiven.add( famId );
					TALE[] mems = families[i].getFamilyMembers();
					for(int j=0;j<mems.length;j++){
						String tid = famId+(j+1);
						reserved.add( tid );
						tid += (orgSuff == null ? "" : " "+orgSuff);
						
						String oldId = mems[j].getId();
						
						if(orgSuff != null){
							mems[j].setStrain(orgSuff);
						}
						if(accession != null){
							mems[j].setAccession(accession);
						}
						
						nameMap.add( new String[]{oldId,tid} );
						if( !mems[j].getId().contains("tempTALE") ){
							tid += " ("+ mems[j].getId() + ")";
						}else if(mems[j].getId().contains("Pseudo") || mems[j].getId().contains("pseudo")){
							tid += " (Pseudo)";
						}
						mems[j].setId( tid );
						
					}
				}else{
					TALE[] mems = families[i].getFamilyMembers();
					for(int j=0;j<mems.length;j++){
						if(mems[j].isNew()){
							String tid = families[i].getFamilyId();
							int k=1;
							String temp = tid+k;
							while(reserved.contains( temp )){
								k++;
								temp = tid+k;
							}
							tid = temp;
							reserved.add( tid );
							tid += (orgSuff == null ? "" : " "+orgSuff);
							
							String oldId = mems[j].getId();
							
							if(orgSuff != null){
								mems[j].setStrain(orgSuff);
							}
							if(accession != null){
								mems[j].setAccession(accession);
							}
							
							nameMap.add( new String[]{oldId,tid} );
							if( !mems[j].getId().contains("tempTALE") ){
								tid += " ("+ mems[j].getId() + ")";
							}else if(mems[j].getId().contains("Pseudo") || mems[j].getId().contains("pseudo")){
								tid += " (Pseudo)";
							}
							mems[j].setId( tid );
						}
					}
				}
			}
			
			String[][] temp = nameMap.toArray( new String[0][0] ); 
			
			String[][] temp2 = new String[temp.length+lastNameMap.length][];
			
			System.arraycopy(lastNameMap, 0, temp2, 0, lastNameMap.length);
			System.arraycopy(temp, 0, temp2, lastNameMap.length, temp.length);
			
			lastNameMap = temp2;
			
			
		}
		
	}
	
	
	private static NumberFormat format = DecimalFormat.getInstance( Locale.US );
	private static NumberFormat formatE = new DecimalFormat("0.##E0");
	
	public ClassAssignmentTool() {
		
	}

	@Override
	public ToolParameterSet getToolParameters() {

		FileParameter builderFile = new FileParameter( "Class builder", "TALE class builder definition", "xml", true );
		builderFile.setExtendedType( TALEFamilyBuilder.class.getName() );
		
		FileParameter fp = new FileParameter("TALE sequences","The sequences of the TALEs (DNA or protein), or \"TALE DNA parts\" or \"TALE Protein parts\" output of \"TALE Analysis\", or RVD sequences.","fasta,fa,fas",true);
		fp.setExtendedType( "fasta/dna" );
		
		SimpleParameter sp = null;
		SimpleParameter ap = null;
		try {
			sp = new SimpleParameter( DataType.STRING, "Strain", "The name of the strain.", false );
			ap = new SimpleParameter( DataType.STRING, "Accession", "The accesion number of the genome (if applicable).", false );
		} catch ( DatatypeNotValidException e ) {
			e.printStackTrace();
		}
		
		return new ToolParameterSet( getShortName(), builderFile, fp, sp, ap );
	}

	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {
		
		progress.setLast( 1.0 );
		progress.setCurrent( 0.0 );
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder(new StringBuffer( ((FileParameter)parameters.getParameterAt( 0 )).getFileContents().getContent()));
		TALE[] newttales = null;
		try{
			newttales = ClassBuilderTool.readProteinTALEs( ((FileParameter)parameters.getParameterAt( 1 )).getFileContents(), protocol );
		}catch(Exception e){
			protocol.appendWarning("TALE sequence could not be read. Please make sure that the input is indeed a FastA containing TALE sequences or TALE parts.\n");
			protocol.appendThrowable(e);
			throw e;
		}
		
		
		String orgSuff = (String)parameters.getParameterAt( 2 ).getValue();
		if(orgSuff != null && orgSuff.trim().length() == 0){
			orgSuff = null;
		}
		String accession = (String)parameters.getParameterAt( 3 ).getValue();
		if(accession != null && accession.trim().length() == 0){
			accession = null;
		}
		
		
		double cut = builder.getCut();
		
		
		builder.setToOld();
		
		TALEFamily[] fams = builder.getFamilies();
		
		LinkedList<Result> topLevel = new LinkedList<Result>();
		LinkedList<Result> ress = new LinkedList<Result>(); 
		
		ArrayList<TALE> notAssigned = new ArrayList<TALE>();
		
		
		boolean[] added = new boolean[fams.length];
		
		//LinkedList<TALE>[] assignments = new LinkedList[fams.length];
		protocol.append( "Assigning TALE...\n" );
		
		SchemaFamilyIdGenerator gen = new SchemaFamilyIdGenerator( orgSuff, accession );
		
		for(int i=0;i<newttales.length;i++){
			protocol.append( newttales[i].getId() + "\n" );
			newttales[i].setIsNew(true);
			
			StringBuffer report = new StringBuffer();

			report.append("New TALE: "+newttales[i].getId()+"\n");

			if(newttales[i].getNumberOfRepeats()>3){
				
				Pair<Integer,Double> idx = builder.getClosestFamilyIndex( newttales[i], null );
				
				/*if(newttales[i].getId().equals("tempTALE9 (Pseudo)") && orgSuff.equals("Xoo PXO142")){
					for(int j=0;j<fams.length;j++){
						if(fams[j].getFamilyId().equals("TalAQ")){
							idx = new Pair<Integer, Double>(j, fams[j].getDistance( newttales[i], null, builder ));
						}
					}
				}*/
				
				
				if(idx.getSecondElement() < cut){
					report.append("Assigned to class "+fams[idx.getFirstElement()].getFamilyId()+".\n");
					double p2 = fams[idx.getFirstElement()].getSignificance(newttales[i], null, null, builder);
					report.append("distance "+format.format(idx.getSecondElement())+", (p="+formatE.format(Math.pow( 10, p2 ))+")\n");
					//newFams[idx.getFirstElement()] = newFams[idx.getFirstElement()].addTALE(newttales[i]);

					/*LinkedList<TALE>[] assignments = new LinkedList[fams.length];
					if(assignments[idx.getFirstElement()] == null){
						assignments[idx.getFirstElement()] = new LinkedList();
					}
					assignments[idx.getFirstElement()].add( newttales[i] );*/
					
					LinkedList<TALE> temp = new LinkedList<TALE>();
					temp.add(newttales[i]);
					
					LinkedList<Pair<Integer,LinkedList<TALE>>> ass = new LinkedList<Pair<Integer,LinkedList<TALE>>>();
					ass.add(new Pair<Integer, LinkedList<TALE>>(idx.getFirstElement(), temp));
					

					builder.addTALEsToFamilies( ass.toArray( new Pair[0] ), new TALE[0], gen);
					
					fams[idx.getFirstElement()] = builder.getFamily(idx.getFirstElement());
					
					newttales[i].setIsNew(false);
					
					added[idx.getFirstElement()] = true;
				}else{
					report.append("Not assigned to any class.\n");
					notAssigned.add( newttales[i] );
					idx = new Pair<Integer, Double>(-1, null);
				}

				report.append("Other families with significant matches:\n");
				boolean found = false;
				for(int j=0;j<fams.length;j++){
					if(j != idx.getFirstElement()){
						double p = fams[j].getSignificance( newttales[i], null, null, builder );
						//System.out.println(fams[j].getFamilyId()+" "+format.format(fams[j].getDistance( newttales[i], null, builder ))+" "+p);
						if(p < Math.log10( 0.001 )){
							double d = fams[j].getDistance( newttales[i], null, builder );
							report.append("Class "+fams[j].getFamilyId()+": distance "+format.format(d)+", (p="+formatE.format(Math.pow( 10, p ))+"),\n");
							found = true;
						}
					}
				}
				if(!found){
					report.append("none\n");
				}
				
				report.append("\nPlease note that numbers within TALE classes are assigned successively and may change when TALEs are included into the official AnnoTALE database later.\n");

				ress.add( new TextResult( "Report for "+newttales[i].getId(), "", new FileRepresentation( "", report.toString() ), "txt", "TALE Class Assignment", null, false ) );

			}else{
				report.append("Not assigned to a class, because it has less than 4 repeats.\n");
				protocol.appendWarning("TALE "+newttales[i].getId()+" not assigned to a class, because it has less than 4 repeats.\n");
				
				ress.add( new TextResult( "Report for "+newttales[i].getId(), "", new FileRepresentation( "", report.toString() ), "txt", "TALE Class Assignment", null, false ) );
			}
			progress.setCurrent( i/(double)newttales.length/2.0 );
		}
		
		topLevel.add( new ResultSetResult( "Reports", "Assignment reports for new TALEs", null, new ResultSet( ress.toArray( new Result[0] ) ) ) );
		
		ress.clear();
		
		/*LinkedList<Pair<Integer,LinkedList<TALE>>> ass = new LinkedList<Pair<Integer,LinkedList<TALE>>>();
		
		for(int i=0;i<assignments.length;i++){
			if(assignments[i] != null){
				ass.add( new Pair<Integer, LinkedList<TALE>>( i, assignments[i] ) );
			}
		}*/
		
		
		
	//	protocol.append( "\nAdding new classes...\n" );
		
		builder.addTALEsToFamilies( new Pair[0], notAssigned.toArray( new TALE[0] ), gen);
		
		progress.setCurrent( 3.0/4.0 );
		
		for(int i=0;i<newttales.length;i++){
			newttales[i].setIsNew(true);
		}
		
		TALEFamily[] newFams = builder.getFamilies();
		
		ComparableElement<Boolean, TALEFamily>[] ces = new ComparableElement[added.length];
		for(int i=0;i<added.length;i++){
			ces[i] = new ComparableElement<Boolean, TALEFamilyBuilder.TALEFamily>( added[i], newFams[i] );
		}
		Arrays.sort( ces );
		for(int i=0;i<added.length;i++){
			added[i] = ces[i].getElement();
			newFams[i] = ces[i].getWeight();
		}
				
		TALgetterDiffSM model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLI.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/talfinder_obg2_hyp_bg.xml" ) ), "model" );
		
		for(int i=0;i<added.length;i++){
			if(added[i]){
				PlotGeneratorResult pgr = new PlotGeneratorResult( "Class tree for "+newFams[i].getFamilyId(), "Plot of the tree of the TALEs in this class", newFams[i], true );
				TextResult fileres = new TextResult( "Class report for "+newFams[i].getFamilyId(), "Report for class "+newFams[i].getFamilyId(), new FileRepresentation( "", newFams[i].toString(model, builder) ), "txt", "TALE Class Assignment", null, false );
				ResultSetResult rsr = new ResultSetResult( "Modified class "+newFams[i].getFamilyId(), "Collection of results for class "+newFams[i].getFamilyId(), null, new ResultSet( new Result[]{fileres,pgr} ) );
				
				ress.add(rsr);
			}
		}
		
		
		if(notAssigned.size()>0){
			protocol.append( "Creating new classes...\n" );
			FamilyResult[] res2 = BuildFamilies.getFamilyResults( newFams, builder.getPVal(), builder, added.length );
			Arrays.sort( res2 );
			protocol.append( "... found "+res2.length+" new classes.\n" );
			for(int i=0;i<res2.length;i++){
				PlotGeneratorResult pgr = new PlotGeneratorResult( "Class tree for "+res2[i].getFamily().getFamilyId(), "Plot of the tree of the TALEs in this class", res2[i].getFamily(), true );
				TextResult fileres = new TextResult( "Class report for "+res2[i].getFamily().getFamilyId(), "Report for class "+res2[i].getFamily().getFamilyId(), new FileRepresentation( "", res2[i].getFamily().toString(model, builder) ), "txt", "TALE Class Assignment", null, false );
				ResultSetResult rsr = new ResultSetResult( "New class "+res2[i].getFamily().getFamilyId(), "Collection of results for class "+res2[i].getFamily().getFamilyId(), null, new ResultSet( new Result[]{fileres,pgr} ) );
				
				ress.add(rsr);
			}
		}
		
		topLevel.add( new ResultSetResult( "Classes", "Modified and new TALE classes", null, new ResultSet( ress.toArray( new Result[0] ) ) ) );
		ress.clear();
		
		protocol.append( "\nProposing TALE names...\n" );
		
		
		StringBuffer summary = new StringBuffer();
		
		
		
		
		LinkedList<ResultSet> props = new LinkedList<ResultSet>();
		
		String[][] map = gen.lastNameMap;
		
		summary.append("The "+map.length+" TALEs have been assigned to the following classes:\n");
		for(int i=0;i<map.length;i++){
			props.add( new ResultSet( new Result[]{
			                                       new CategoricalResult( "Old name", "", map[i][0] ),
			                                       new CategoricalResult( "New name", "", map[i][1] )
			} ) );
			String clazz = map[i][1];
			String temp = clazz.replaceAll("\\s.*$", "");
			clazz = temp.replaceAll("[0-9]+$", "");
			int num = Integer.parseInt(temp.substring(clazz.length()));
			
			summary.append("- "+map[i][0]+" has been assigned to "+(num == 1 ? "new" : "existing")+" class "+clazz+" with name "+map[i][1]+"\n");
		}
		
		String oss = orgSuff == null ? "" : " ("+orgSuff+")";
		
		summary.append("\nPlease note that numbers within TALE classes are assigned successively and may change when TALEs are included into the official AnnoTALE database later.\n");
		
		summary.append("\n\nA detailed assignment report for each individual TALE can be found under \"Reports\".\n\n"
				+ "Class reports and class trees for all classes that have been created or \n"
				+ "modified due to class assignment of these TALEs are available under \"Classes\".\n\n"
				+ "The systematic names proposed for these TALEs based on the class assignment are \n"
				+ "available as \""+"TALE names"+oss+"\".\n"
				+ "This file can also be used for the \"Rename TALEs in File\" module of AnnoTALE.\n\n"
				+ "A class builder that may be used in successive runs of the \"TALE Class Assignment\" \n"
				+ "module is provided as \""+"Augmented class builder"+oss+"\".\n\n"
				+ "Renamed TALE DNA and protein sequences are provided under \""+"Renamed TALE protein sequences"+oss+"\".");
		
		
		topLevel.addFirst(new TextResult("Summary of class assignment", "A summary of the class assignment", new FileRepresentation("", summary.toString()), "txt", "TALE Class Builder", TALEFamilyBuilder.class.getName(), true));
		
		
		
		
		
		ListResult lr = new ListResult( "TALE names"+oss, "Proposed TALE names based on class assignment", null, props.toArray( new ResultSet[0] ) );
		
	//	topLevel.add( new TextResult( "TALE names", "Proposed TALE names based on class assignment", new FileRepresentation( "", buf.toString() ), "tsv", "TALE Class Builder", null ) );
		
		topLevel.add( lr );
		
		topLevel.add( new TextResult( "Augmented class builder"+oss, "TALE class builder definition", new FileRepresentation( "", builder.toXML().toString() ), "xml", "TALE Class Builder", TALEFamilyBuilder.class.getName(), true ) );
		
		
		
		
		StringBuffer dna = new StringBuffer();
		StringBuffer prot = new StringBuffer();
		protocol.append( "Writing FastA outputs.\n" );
		for(int i=0;i<newttales.length;i++){
			TALE dnaTALE = newttales[i].getDnaOriginal();
			if(dnaTALE != null){
				dna.append( ">"+dnaTALE.getId()+"\n" );
				dna.append( dnaTALE.getStart() );
				for(int j=0;j<dnaTALE.getNumberOfRepeats();j++){
					dna.append( dnaTALE.getRepeat( j ).getRepeat() );
				}
				dna.append( newttales[i].getDnaOriginal().getEnd()+"\n" );
			}
			
			prot.append( ">"+newttales[i].getId()+"\n" );
			prot.append( newttales[i].getStart() );
			for(int j=0;j<newttales[i].getNumberOfRepeats();j++){
				prot.append( newttales[i].getRepeat( j ).getRepeat() );
			}
			prot.append( newttales[i].getEnd()+"\n" );
		}
		
		TextResult fres4 = new TextResult( "Renamed TALE protein sequences"+oss, "The protein sequences of the assigned TALEs using the proposed names", new FileRepresentation( "", prot.toString() ), "fasta", "Class Assignment", "fasta/as", true );
		
		ResultSetResult seqs = null;
		
		if(dna.length() > 0){
			TextResult fres3 = new TextResult( "Renamed TALE DNA sequences"+oss, "The DNA sequences of the assigned TALEs using the proposed names", new FileRepresentation( "", dna.toString() ), "fasta", "Class Assignment", "fasta/dna", true );
			seqs = new ResultSetResult( "Renamed TALE sequences"+oss, "The sequences of the assigned TALEs using the proposed names", null, new ResultSet( new Result[]{fres3,fres4} ) );
		}else{
			seqs = new ResultSetResult( "Renamed TALE sequences"+oss, "The sequences of the assigned TALEs using the proposed names", null, new ResultSet( new Result[]{fres4} ) );
		}
		
		topLevel.add( seqs );
		
		progress.setCurrent( 1.0 );
		
		ResultSet set = new ResultSet( topLevel.toArray( new Result[0] ) );
		return new ToolResult("Result of "+getToolName()+oss, getToolName()+" on "+((FileParameter)parameters.getParameterAt( 1 )).getFileContents().getFilename(), null, set, parameters, getToolName(), new Date(System.currentTimeMillis()));
	}

	@Override
	public String getToolName() {
		return "TALE Class Assignment";
	}

	@Override
	public String getShortName() {
		return "assign";
	}

	@Override
	public String getDescription() {
		return "Assigns TALEs to a TALE class and proposes names";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( ClassAssignmentTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/ClassAssignmentTool.txt" ) ).toString();
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

}
