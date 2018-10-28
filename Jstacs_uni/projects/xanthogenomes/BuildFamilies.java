package projects.xanthogenomes;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Locale;

import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.StringAlignment;
import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.clustering.hierachical.Hclust.Linkage;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.SeqLogoPlotter;
import de.jstacs.utils.graphics.GraphicsAdaptor;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory;
import de.jstacs.utils.graphics.GraphicsAdaptorFactory.OutputFormat;
import projects.tals.ScanForTBSCLI;
import projects.tals.TALgetterDiffSM;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;
import projects.xanthogenomes.Tools.Translator;
import projects.xanthogenomes.alignmentCosts.RVDCosts;


public class BuildFamilies {

	private static double extraGapOpening = 1.0;
	private static double extraGapExtension = 0.1;
	private static Linkage linkage = Linkage.AVERAGE;
	private static AlignmentType at = AlignmentType.SEMI_GLOBAL;
	
	private static NumberFormat format = DecimalFormat.getInstance( Locale.US );
	private static NumberFormat formatE = new DecimalFormat("0.##E0");
	
	public static class FamilyResult implements Comparable<FamilyResult> {
		
		private TALEFamily family;
		private TALEFamily[] relatedFams;
		private TALE[][] relatedTALEs;
		private int[][] myTALEs;
		private StringAlignment[][] alignments;
		private double[][] pvals;
		
		private ClusterTree<TALEFamily> familyTree;
		
		/**
		 * @param family
		 * @param relatedFams
		 * @param relatedTALEs
		 * @param relationScores
		 */
		public FamilyResult( TALEFamily family, TALEFamily[] relatedFams, TALE[][] relatedTALEs, StringAlignment[][] alignments, double[][] pvals, int[][] myTALEs ) {
			this.family = family;
			this.relatedFams = relatedFams;
			this.relatedTALEs = relatedTALEs;
			this.alignments = alignments;
			this.myTALEs = myTALEs;
			this.pvals = pvals;
		}
		
		public int compareTo(FamilyResult fr){
			return family.getFamilyId().compareTo( fr.family.getFamilyId() );
		}
		
		public void setSimilarityGroup(ClusterTree<TALEFamily>[] famTrees){
			for(int i=0;i<famTrees.length;i++){
				TALEFamily[] fams = famTrees[i].getClusterElements();
				for(int j=0;j<fams.length;j++){
					if(fams[j].getFamilyId().equals( this.family.getFamilyId() )){
						this.familyTree = famTrees[i];
						return;
					}
				}
			}
		}
		
		public TALEFamily getFamily() {
			return family;
		}
		
		public TALEFamily[] getRelatedFams() {
			return relatedFams;
		}
		
		public TALE[][] getRelatedTALEs() {
			return relatedTALEs;
		}
		
		public StringAlignment[][] getAlignments() {
			return alignments;
		}
		
		public String toString(TALgetterDiffSM model, TALEFamilyBuilder builder){
			StringBuffer sb = new StringBuffer();
			sb.append( family.toString(model,builder)+"\n" );
			if(relatedFams.length == 0){
				sb.append( "No related classes.\n" );
			}else{
				TALE[] members = family.getFamilyMembers();
				sb.append( "Related classes with significant matches:\n" );
				for(int i=0;i<relatedFams.length;i++){
					sb.append( "Class "+relatedFams[i].getFamilyId()+" with members\n" );
					for(int j=0;j<relatedTALEs[i].length;j++){
						sb.append( relatedTALEs[i][j].getId()+" related to "+members[ myTALEs[i][j] ].getId()+" with score "+format.format( alignments[i][j].getCost() )+" (p="+formatE.format( Math.pow(10,pvals[i][j]) )+")\n");
					}
					sb.append( "\n" );
					sb.append( "Alignments:\n" );
					for(int j=0;j<alignments[i].length;j++){
						sb.append(relatedTALEs[i][j].getId()+" vs. "+members[ myTALEs[i][j] ].getId()+":\n");
						sb.append( TALEAligner.alignmentToString( alignments[i][j], format )+"\n\n" );
					}
				}
			}
			sb.append( "\n" );
			if(this.familyTree != null){
				if(this.familyTree.getNumberOfElements() == 1){
					sb.append( "This class belongs to its own similarity group\n" );
				}else{
					sb.append( "This class belongs to a similarity group together with classes " );
					TALEFamily[] group = this.familyTree.getClusterElements();
					for(int j=0;j<group.length;j++){
						if(!group[j].getFamilyId().equals( this.family.getFamilyId() )){
							sb.append( group[j].getFamilyId()+", " );
						}
					}
					sb.delete( sb.length()-2, sb.length() );
					sb.append( "\n" );
				}
			}
			
			
			return sb.toString();
		}
		
		
		
	}
	
	
	public static FamilyResult[] getFamilyResults(TALEFamily[] fams, double pval, TALEFamilyBuilder builder, int offset){
		
		Costs costs = builder.getCosts();
		
		AlignmentPValues pv = null;
		
		if(costs instanceof AffineCosts && ((AffineCosts)costs).getInternalCosts() instanceof RVDCosts){
			
			pv = new AlignmentPValues( builder.getAllTALEs(), (RVDCosts)((AffineCosts)costs).getInternalCosts() );
		}
		
		
		
		FamilyResult[] ress = new FamilyResult[fams.length-offset];
		
		for(int i=offset;i<fams.length;i++){
			TALE[] members = fams[i].getFamilyMembers();
			LinkedList<TALEFamily> relFams = new LinkedList<TALEFamilyBuilder.TALEFamily>();
			LinkedList<TALE[]> relsTales = new LinkedList<TALE[]>();
			LinkedList<StringAlignment[]> relsAls = new LinkedList<StringAlignment[]>();
			LinkedList<int[]> relsMys = new LinkedList<int[]>();
			LinkedList<double[]> relPs = new LinkedList<double[]>();
			for(int j=0;j<fams.length;j++){
				if(i != j){
					
					LinkedList<TALE> relTales = new LinkedList<TALE>();
					LinkedList<StringAlignment> relAls = new LinkedList<StringAlignment>();
					IntList relMy = new IntList();
					DoubleList relP = new DoubleList();
					TALE[] members2 = fams[j].getFamilyMembers();
					for(int k=0;k<members.length;k++){
						for(int l=0;l<members2.length;l++){
							
							StringAlignment al = TALEAligner.align( members[k], members2[l], costs, at, extraGapOpening, extraGapExtension );
							if(pv != null){
								double p = pv.getLog10PValue( members[k], members2[l], al.getCost(), extraGapOpening, extraGapExtension );
								if(p<Math.log10( pval )){
									relTales.add( members2[l] );
									relAls.add( al );
									relMy.add( k );
									relP.add( p );
								}
							}
							
						}
					}
					
					if(relTales.size() > 0){
						relFams.add( fams[j] );
						relsTales.add( relTales.toArray( new TALE[0] ) );
						relsAls.add( relAls.toArray( new StringAlignment[0] ) );
						relsMys.add( relMy.toArray() );
						relPs.add( relP.toArray() );
					}
				}
								
			}
			
			ress[i-offset] = new FamilyResult( fams[i], relFams.toArray( new TALEFamily[0] ), relsTales.toArray( new TALE[0][0] ), relsAls.toArray( new StringAlignment[0][0] ), relPs.toArray( new double[0][0] ), relsMys.toArray( new int[0][0] ) );
			
		}
		
		return ress;
	}
	
	
	public static Pair<TALEFamilyBuilder,FamilyResult[]> build(TALE[] ttales, double cut, double pval) throws NonParsableException, IOException{
		RVDCosts rvdCosts = new RVDCosts( 1.0, 0.2, 0.8, 0.0 );//TODO
		Costs costs = new AffineCosts( 5.0, 5.0, rvdCosts );
		
		//TALgetterDiffSM model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLI.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/talfinder_obg2_hyp_bg.xml" ) ), "model" );

		//SpecificityCosts specCosts = new SpecificityCosts( 1.0, model );//TODO
		//Costs costs = new AffineCosts( 5.0, specCosts );
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder( ttales, costs, linkage, at, extraGapOpening, extraGapExtension, cut, pval );
		
		//AlignmentPValues pv = null;//new AlignmentPValues( ttales, rvdCosts );//TODO
		
		TALEFamily[] fams = builder.getFamilies();
		
		FamilyResult[] ress = getFamilyResults( fams, pval, builder, 0 );
		
		return new Pair<TALEFamilyBuilder,FamilyResult[]>(builder,ress);
		
	}
	
	/*public static Pair<ClusterTree<TALEFamily>[],ClusterTree<TALEFamily>> groupBySpecificity(TALgetterDiffSM model, TALEFamilyBuilder builder, double cut){
		SpecificityCosts specCosts = new SpecificityCosts( 1.0, model );
		Costs specCosts2 = new AffineCosts( 5.0, specCosts );
		
		TALEFamily[] specFams = builder.rebuild( linkage, specCosts2, at, extraGapOpening, extraGapExtension );
		
		ClusterTree<TALEFamily> famTree = TALEFamilyBuilder.clusterFamilies( specFams, linkage, specCosts2, at, extraGapOpening, extraGapExtension );
	
		ClusterTree<TALEFamily>[] subFamTrees = Hclust.cutTree( cut, famTree );
		
		return new Pair<ClusterTree<TALEFamily>[],ClusterTree<TALEFamily>>(subFamTrees,famTree);
	}*/
	
	
	public static void rename(TALE[] ttales) throws IOException{
		HashMap<String, String> nameMap = new HashMap<String, String>();
		
		BufferedReader read = new BufferedReader( new FileReader( "/Users/dev/Desktop/TAL-Chips/Genomes/proteins/tal_names.csv" ) );//TODO
		
		String str = null;
		while( (str = read.readLine()) != null ){
			String[] parts = str.split( "\t" );
			nameMap.put( parts[0].trim(), parts[2].trim()+" ("+parts[1].trim()+")" );
		}
		read.close();
		
		for(int i=0;i<ttales.length;i++){
			if(nameMap.get( ttales[i].getId() ) == null){
				System.out.println(ttales[i].getId());
			}
			ttales[i].setId( nameMap.get( ttales[i].getId() ) );
		}
	}
	
	
	public static void plotGroups(GraphicsAdaptor adaptor, ClusterTree<TALEFamily> tree, String nameBase) throws IOException{
		FamilyGroupPlotter plotter = new FamilyGroupPlotter( 100 );
		
		Graphics2D dummy = adaptor.getGraphics( 10, 10 );
		
		int[] dim = plotter.getDimension( dummy, tree );
		
		Graphics2D graphics = adaptor.getGraphics( dim[0], dim[1] );
		
		graphics.setColor( Color.white );
		
		graphics.fillRect( 0, 0, dim[0], dim[1] );
		
		graphics.setColor( Color.black );
		
		graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		graphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		
		plotter.plot( graphics, tree );
		
		adaptor.generateOutput( nameBase+"."+adaptor.getGraphicsExtension() );
	}
	
	public static void main( String[] args ) throws Exception {
		GraphicsAdaptor adaptor = GraphicsAdaptorFactory.getAdaptor( OutputFormat.PDF );
		
		//AlphabetContainer protAlph = Tools.Translator.DEFAULT.getProteinAlphabet();
				
		TALE[] tales = TALE.readTALEs( args[0], args[1], args[2] );
		
		String outpath = args[3];
		
		TALE[] ttales = TALE.translateTALEs( tales, Translator.DEFAULT );
		
		rename( ttales );
		
		double cut = 6.0;
		double pval = 0.01;
	//	double cutFam = 6.0;
		
		Pair<TALEFamilyBuilder,FamilyResult[]> res = build( ttales, cut, pval );
		
		FamilyResult[] famRes = res.getSecondElement();
		
		TALgetterDiffSM model = (TALgetterDiffSM)XMLParser.extractObjectForTags( FileManager.readInputStream( ScanForTBSCLI.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/talfinder_obg2_hyp_bg.xml" ) ), "model" );
		
		//Pair<ClusterTree<TALEFamily>[],ClusterTree<TALEFamily>> groups = groupBySpecificity( model, res.getFirstElement(), cutFam );
		
		ClusterTree<TALEFamily> famTree = res.getFirstElement().clusterFamilies();
		
		plotGroups( adaptor, famTree, outpath+"/family_tree" );
		
		for(int i=0;i<famRes.length;i++){
		//	famRes[i].setSimilarityGroup( groups.getFirstElement() );
			TALEFamily fam = famRes[i].getFamily();
			fam.plotFamilyToFile(outpath+"/family_"+(i+1), adaptor );
			PrintWriter wr = new PrintWriter( outpath+"/family_"+(i+1)+"_report.txt" );
			wr.println( famRes[i].toString(model, res.getFirstElement()) );
			wr.close();
			double[][] specs = fam.getSpecificityProfile( model );
			int w = SeqLogoPlotter.getWidth( 300, specs );
			Graphics2D graph = adaptor.getGraphics( w, 300 );
			String[] labels = new String[specs.length];
			for(int j=0;j<labels.length;j++){
				labels[j] = j+"";
			}
			SeqLogoPlotter.plotLogo( graph, 300, specs, labels, "Position", "bits" );
			adaptor.generateOutput( outpath+"/family_"+(i+1)+"_specificity."+adaptor.getGraphicsExtension() );
		}
		
		

	}

}
