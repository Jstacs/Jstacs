package projects.xanthogenomes.tools;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.io.FileManager;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import projects.talen.InfixMatchFinder;
import projects.talen.MatchFinder;
import projects.talen.MatchFinder.Match;
import projects.tals.ScanForTBSCLI;
import projects.tals.TALgetterDiffSM;
import projects.tals.TBSScanner;
import projects.xanthogenomes.RVDAlphabetContainer;
import projects.xanthogenomes.TALE;
import projects.xanthogenomes.TALEFamilyBuilder;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;


public class PredictAndIntersectTargetsTool implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		FileParameter genome = new FileParameter( "Input sequences", "Sequences, e.g., promoters, to scan for TALE target sites", "fasta,fa,fas", true );
		
		
		FileParameter builderFile = new FileParameter( "Class builder", "TALE class builder definition", "xml", true );
		builderFile.setExtendedType( TALEFamilyBuilder.class.getName() );
		
		FileParameter tales = new FileParameter( "TALE sequences", "TALE sequences, either as complete DNA or AS sequences (e.g., output of TALE Prediction) or as RVD sequences.", "fasta,fa,fas",true);
		
		
		SelectionParameter selPar = null;
		try{
			selPar = new SelectionParameter( DataType.PARAMETERSET, new String[]{"TALEs in FastA","TALEs in class builder"}, new ParameterSet[]{new SimpleParameterSet( tales ), new SimpleParameterSet( builderFile )}, null, "Predictions for", "Predict and intersect targets for all TALEs in a given FastA input file or for all TALEs in all classes defined by a class builder.", true );
		}catch(Exception e){
			e.printStackTrace( );
			return null;
		}
		
		
		return new ToolParameterSet( getShortName(), genome,selPar);
	}

	@Override
	public ToolResult run( ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads ) throws Exception {

		progress.setLast( 1.0 );
		progress.setCurrent( 0.0 );
		
		FileParameter fp = (FileParameter)parameters.getParameterAt( 0 );
		FileRepresentation fr = fp.getFileContents();
		
		protocol.append( "Reading input data...\n" );
		
		DataSet ds = new DataSet( new AlphabetContainer(new DiscreteAlphabet( true, "A", "C", "G", "T", "N", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V" )), new SparseStringExtractor( new StringReader( fr.getContent() ), '>', "",new SimpleSequenceAnnotationParser() ) ); 
		
		
		
		Pair<int[][],DataSet> pair = TBSScanner.preprocess( ds );

		int[][] offsets = pair.getFirstElement();
		ds = pair.getSecondElement();

		progress.setCurrent( 0.1 );
		protocol.append( "...finished.\n\n" );
		
		int cap = 500;
		
		SelectionParameter selPar = (SelectionParameter)parameters.getParameterAt( 1 );
		
		
		
		StringBuffer[] sbs = null;
		String[] classNames = null;
		
		protocol.append( "Collecting TALE RVD sequences...\n" );
		
		TALE[][] allTales = null;
		
		if(selPar.getSelected() == 0){

			StringBuffer sb = null;
			try{
				TALE[] tales = ClassBuilderTool.readProteinTALEs( ( (FileParameter)((ParameterSet)selPar.getValue()).getParameterAt( 0 ) ).getFileContents() , protocol );
				
				tales = filterTALEsByLength(tales, protocol);
				sb = new StringBuffer();
				for(int i=0;i<tales.length;i++){
					sb.append( ">"+tales[i].getId()+"\n" );
					Sequence rvds = tales[i].getRvdSequence();
					sb.append( rvds.toString( "-", 0, rvds.getLength() )+"\n" );
				}

				allTales = new TALE[][]{tales};
			}catch(Exception e){
				sb = new StringBuffer();
				sb.append( ( (FileParameter)((ParameterSet)selPar.getValue()).getParameterAt( 0 ) ).getFileContents().getContent() );
			}

			sbs = new StringBuffer[]{sb};
		}else{
			
			TALEFamilyBuilder builder = new TALEFamilyBuilder(new StringBuffer( ((FileParameter)((ParameterSet)selPar.getValue()).getParameterAt( 0 )).getFileContents().getContent()));

			TALEFamily[] fams = builder.getFamilies();
			Arrays.sort( fams );
			classNames = new String[fams.length];
			
			sbs = new StringBuffer[fams.length];
			allTales = new TALE[fams.length][];
			
			for(int j=0;j<fams.length;j++){
				classNames[j] = fams[j].getFamilyId();
				TALE[] tales = fams[j].getFamilyMembers();
				tales = filterTALEsByLength(tales, protocol);
				allTales[j] = tales;
				StringBuffer sb = new StringBuffer();
				for(int i=0;i<tales.length;i++){
					sb.append( ">"+tales[i].getId()+"\n" );
					Sequence rvds = tales[i].getRvdSequence();
					sb.append( rvds.toString( "-", 0, rvds.getLength() )+"\n" );
				}
				sbs[j] = sb;
			}
			
		}
		
		

		
		TALgetterDiffSM model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLI.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/talfinder_obg2_hyp_bg.xml" ) ), "model" );
		DiscreteAlphabet modelAlph = (DiscreteAlphabet)model.getRVDAlphabet().getAlphabetAt( 0 );
		HashSet<String> newSyms = new HashSet<String>();
		
		double totNum = 0.0;
		
		
		for(int t=0;t<sbs.length;t++){

			DataSet tals = new DataSet(RVDAlphabetContainer.SINGLETON,new SparseStringExtractor( new StringReader( sbs[t].toString() ), '>', "", new SimpleSequenceAnnotationParser() ),"-");
			for(int i=0;i<tals.getNumberOfElements();i++){
				totNum++;
				Sequence rvds = tals.getElementAt( i );
				for(int j=0;j<rvds.getLength();j++){
					String sym = rvds.toString( j, j+1 );
					if(!modelAlph.isSymbol( sym ) && !newSyms.contains( sym )){
						newSyms.add(sym);
					}
				}
			}
		}
		
		String[] ns = newSyms.toArray( new String[0] );
		double[][] specs = new double[ns.length][];
		for(int i=0;i<ns.length;i++){
			specs[i] = new double[]{0.25,0.25,0.25,0.25};
		}
		model.addAndSet( ns, specs, null );
		
		
		model.fix();

		progress.setCurrent( 0.2 );
		protocol.append( "...finished.\n\n" );
		
		double fac = 1.0/totNum*0.8;
		totNum = 0;
		
		Result[] res = new Result[sbs.length];
		
		ResultSet set = null;
		
		protocol.append( "Predicting targets for\n" );
		for(int s=0;s<sbs.length;s++){

			DataSet tals = new DataSet(model.getRVDAlphabet(),new SparseStringExtractor( new StringReader( sbs[s].toString() ), '>', "", new SimpleSequenceAnnotationParser() ),"-");


			//Sequence tal = Sequence.create( model.getRVDAlphabet(), "NI-NG-NN-NG-NK-NG-NI-NN-NI-NN-NI-NN-NS-NG-NS-NN-NI-N*-NS-NG", "-" );

			Match[][] allMatches = new Match[tals.getNumberOfElements()][];

			String[] talNames = new String[tals.getNumberOfElements()];

			int off = tals.getNumberOfElements() > 1 ? 1 : 0;

			ListResult[] lires = new ListResult[tals.getNumberOfElements()+off];

			for(int t=0;t<tals.getNumberOfElements();t++,totNum++){

				Sequence tal = tals.getElementAt( t );
				talNames[t] = (String)tal.getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue();
				protocol.append( talNames[t]+"\n" );

				//System.out.println("predicting for "+talNames[t]+" "+tal);
				
				//double bestTotal = model.getBestPossibleScore( tal, null );

				ComparableElement<Match, Double>[] list = predict(model, tal, cap, ds, offsets);
				
				if(allTales != null){
					if(allTales[s][t].containsAberrantRepeat()){
					//	System.out.println("is aberrant");
						Sequence tal2 = removeAbberant(tal,allTales[s][t]);
						if(tal2.getLength()>3){
							//	System.out.println("w/ removed: "+tal2);
							ComparableElement<Match, Double>[] list2 = predict(model, tal2, cap, ds, offsets);

							ComparableElement<Match, Double>[] list3 = join(tal,tal2,list,list2,cap);

							list = list3;
						}
					}
				}
				
				
				
				allMatches[t] = new Match[list.length];
				ResultSet[] re = new ResultSet[list.length];
				for(int i=0;i<list.length;i++){
					Match m = list[list.length-1-i].getElement();
					allMatches[t][i] = m;
					double score = list[list.length-1-i].getWeight();
					int seqIdx = m.getSeqIdx();
					int pos = m.getSeqPos()+offsets[0][seqIdx];
					String id = ds.getElementAt( seqIdx ).getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue().toString().trim();

					Sequence currTal = m.getTal() == null ? tal : m.getTal();
					
					//try{
					Sequence ts = ds.getElementAt( seqIdx ).getSubSequence( m.getSeqPos(), currTal.getLength()+1 );
					
					
					
					re[i] = new ResultSet(new Result[]{
					                                   new CategoricalResult( "Sequence ID", "", id ),
					                                   new NumericalResult( "Position", "", pos ),
					                                   new NumericalResult( "Score", "", score ),
					                                   new CategoricalResult( "Site", "", ts.toString() ),
					                                   new CategoricalResult( "Match string", "", model.getMatchString( currTal, ts )+(m.getTal() == null ? "" : "-a") )
					});
				
					/*}catch(IllegalArgumentException e){
						System.out.println(m.getSeqIdx()+" "+m.getSeqPos()+" "+ds.getElementAt(seqIdx).getLength()+" "+m.getTal()+" "+currTal+" "+tal);
						throw e;
					}*/
				}

				lires[t+off] = new ListResult( "Predictions for "+talNames[t] , "TALgetter predictions of target sites for "+talNames[t], null, re );
				progress.setCurrent( 0.2 + fac*totNum );
			}

			if(off > 0){
				ListResult intersection = intersect(ds, offsets,talNames, allMatches);
				lires[0] = intersection;
			}
			//System.out.println(Arrays.toString( lires ));

			
			if(sbs.length == 1){
				set =  new ResultSet( lires );
			}else{
				res[s] = new ResultSetResult( "Predictions for class "+classNames[s], "Predicted target sites for all TALEs in class "+classNames[s], null, new ResultSet( lires ) );
			}

		}
		
		if(set == null){
			set = new ResultSet( res );
		}
		
		progress.setCurrent( 1.0 );
		
		return new ToolResult("Result of "+getToolName(), getToolName()+" on \""+fr.getFilename()+"\"", null, set, parameters, getToolName(), new Date(System.currentTimeMillis()) );
	}
	
	
	private TALE[] filterTALEsByLength(TALE[] tales, Protocol protocol){
		ArrayList<TALE> temp = new ArrayList<TALE>();
		for(int i=0;i<tales.length;i++){
			if(tales[i].getNumberOfRepeats()>3){
				temp.add(tales[i]);
			}else{
				protocol.appendWarning("TALE "+tales[i].getId()+" ignored as it has less than 4 repeats.\n");
			}
		}
		return temp.toArray(new TALE[0]);
	}
	
	private static ComparableElement<Match, Double>[] join(Sequence tal, Sequence tal2, ComparableElement<Match, Double>[] list, ComparableElement<Match, Double>[] list2, int cap) {
		double correct = (tal.getLength()-tal2.getLength())*Math.log(0.25);
		
		
		ComparableElement<Match, Double>[] list3 = new ComparableElement[cap]; 
		
		
		for(int k=list3.length-1,i=list.length-1,j=list2.length-1;k>=0;k--){
			if(i==0){
				list3[k] = list2[j];
				list3[k].getElement().setTal(tal2);
				j--;
			}else if(j==0){
				list3[k] = list[i];
				i--;
			}else{
				double d1 = list[i].getWeight();
				double d2 = list2[j].getWeight();
				d2 += correct;
				if(d1 > d2){
					list3[k] = list[i];
					i--;
				}else{
					list3[k] = new ComparableElement<MatchFinder.Match, Double>(list2[j].getElement(), d2);
					list3[k].getElement().setTal(tal2);
					j--;
				}
			}
			
		}
		
		return list3;
		
	}

	private Sequence removeAbberant(Sequence tal, TALE tale) throws IllegalArgumentException, WrongAlphabetException {
		String temp = tal.toString("-", 0, tal.getLength() );
		String[] parts = temp.split("-");
		StringBuffer tal2 = new StringBuffer();
		for(int i=0;i<tale.getNumberOfRepeats();i++){
			if(tale.getRepeat(i).getType() == TALE.Type.UNKNOWN || tale.getRepeat(i).getType() == TALE.Type.NORMAL || i==tale.getNumberOfRepeats()-1){
				tal2.append(parts[i]);
				tal2.append("-");
			}
		}
		if(tal2.length() > 0){
			tal2.delete(tal2.length()-1, tal2.length());
		}
		
		Sequence res = Sequence.create(tal.getAlphabetContainer(), tal2.toString(), "-");
		return res;
	}

	private static ComparableElement<Match, Double>[] predict(TALgetterDiffSM model, Sequence tal, int cap, DataSet ds, int[][] offsets){
		InfixMatchFinder singleFind = new InfixMatchFinder( null, Math.min( 8, tal.getLength() ), model );


		ComparableElement<Match, Double>[] list = new ComparableElement[0];

		double bestRelScore = model.getBestPossibleScore( tal, null )/(tal.getLength()+1);

		double rat = 0.5;

		while(list.length < cap && rat > 0.1){

			//	System.out.println(model.getBestPossibleScore( tal, null ));
			double singleThresh =  bestRelScore+ Math.log( rat );
			//	System.out.println(singleThresh);

			singleFind.getPreps( tal, singleThresh*(tal.getLength()+1) );

			singleFind.setDataSet( ds );

			list = singleFind.getScoresAbove( tal, singleThresh*(tal.getLength()+1), cap, true, false ).getSortedList();

			rat /= 1.5;

		}

		
		
		return list;
	}
	
	
	private static class Inter{
		
		String id;
		IntList[] ranks;
		IntList[] pos;
		
		public Inter(String id, int len){
			this.id = id;
			this.ranks = new IntList[len];
			this.pos = new IntList[len];
		}
		
		public void add(int tal, int rank, int position){
			if(ranks[tal] == null){
				ranks[tal] = new IntList();
			}
			if(pos[tal] == null){
				pos[tal] = new IntList();
			}
			ranks[tal].add( rank );
			pos[tal].add( position );
		}
		
		public int getLength(){
			return ranks.length;
		}
		
		public int getNum(){
			int num = 0;
			for(int i=0;i<ranks.length;i++){
				if(ranks[i] != null){
					num++;
				}
			}
			return num;
		}
		
		public double getRankProduct(){
			double rp = 1;
			double n = 0;
			for(int i=0;i<ranks.length;i++){
				if(ranks[i] != null){
					double temp = 1.0;
					for(int j=0;j<ranks[i].length();j++){
						temp *= ranks[i].get( j );
					}
					temp = Math.pow( temp, 1.0/ranks[i].length() );
					rp *= temp;
					n++;
				}
			}
			return Math.pow( rp, 1.0/n );
		}

		public String getString( int j ) {
			if(ranks[j] == null){
				return "";
			}
			StringBuffer sb = new StringBuffer();
			for(int i=0;i<ranks[j].length();i++){
				if(i>0){
					sb.append( "; " );
				}
				sb.append( "("+ranks[j].get( i )+","+pos[j].get( i )+")" );
			}
			return sb.toString();
		}
		
	}
	

	private ListResult intersect( DataSet ds, int[][] offsets, String[] talNames, Match[][] allMatches ) {
		HashMap<String, LinkedList<int[]>> map = new HashMap<String, LinkedList<int[]>>();
		
		for(int i=0;i<allMatches.length;i++){
			for(int j=0;j<allMatches[i].length;j++){
				Match m = allMatches[i][j];
				
				int seqIdx = m.getSeqIdx();
				int pos = m.getSeqPos()+offsets[0][seqIdx];
				String id = (String)ds.getElementAt( seqIdx ).getSequenceAnnotationByType( "unparsed comment line", 0 ).getResultForName( "unparsed comment" ).getValue().toString().trim();
							
				if(!map.containsKey( id )){
					LinkedList<int[]> li = new LinkedList<int[]>();
					map.put( id, li );
				}
				map.get( id ).add( new int[]{i,pos,j+1} );
			}
		}
		
		Iterator<String> matches = map.keySet().iterator();
		Inter[] res = new Inter[map.size()];
		//String[] keys = new String[map.size()];
		int k = 0;
		while(matches.hasNext()){
			String m = matches.next();
			//keys[k] = m;
			Inter inter = new Inter( m, allMatches.length );
			LinkedList<int[]> li = map.get( m );
			
			for(int i=0;i<li.size();i++){
				int[] t = li.get( i );
				inter.add( t[0], t[2], t[1] );
			}
			res[k] = inter;
			k++;
		}
		
		
		Arrays.sort( res, new Comparator<Inter>() {

			@Override
			public int compare( Inter o1, Inter o2 ) {
				int c = -(new Integer(o1.getNum()).compareTo( new Integer(o2.getNum() ) ) );
				if(c==0){
					Double rp1 = o1.getRankProduct();
					Double rp2 = o2.getRankProduct();
					return rp1.compareTo( rp2 );
				}else{
					return c;
				}
			}
			
		} );
		
		
		LinkedList<ResultSet> rsl = new LinkedList<ResultSet>();
		
		for(int i=0;i<res.length;i++){
			int num = res[i].getNum();
			String id = res[i].id;
			int len = res[i].getLength();
			Result[] rs = new Result[len+2];
			rs[0] = new CategoricalResult( "Sequence ID", "", id );
			rs[1] = new NumericalResult( "Intersection size", "", num );
			for(int j=0;j<len;j++){
				rs[j+2] = new CategoricalResult( talNames[j], "", res[i].getString(j) );
			}
			rsl.add( new ResultSet( rs ) );
		}
		
		
		return new ListResult( "Overlapping target sites", "Overlapping target sites between TALEs according to TALgetter predictions", null, rsl.toArray( new ResultSet[0] ) );
	}

	@Override
	public String getToolName() {
		return "Predict and Intersect Targets";
	}

	@Override
	public String getShortName() {
		return "targets";
	}

	@Override
	public String getDescription() {
		return "Predicts target sites using TALgetter and intersects targets between different TALEs";
	}

	@Override
	public String getHelpText() {
		try {
			return FileManager.readInputStream( TALEPredictionTool.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/tools/PredictAndIntersectTargetsTool.txt" ) ).toString();
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
