package projects.xanthogenomes;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Locale;

import projects.tals.TALgetterDiffSM;
import projects.xanthogenomes.alignmentCosts.RVDCosts;
import projects.xanthogenomes.tools.ClassAssignmentTool;
import de.jstacs.Storable;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.StringAlignment;
import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.Costs;
import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.clustering.hierachical.Hclust;
import de.jstacs.clustering.hierachical.Hclust.Linkage;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.PlotGeneratorResult.PlotGenerator;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.graphics.GraphicsAdaptor;


public class TALEFamilyBuilder implements Storable {

	private static NumberFormat format = DecimalFormat.getInstance( Locale.US );
	private static NumberFormat formatE = new DecimalFormat("0.##E0");
	
	enum FamilyDistance{
		MIN,
		MAX,
		MEAN
	}
	
	private class DefaultFamilyIdGenerator implements FamilyIdGenerator{

		@Override
		public void setFamilyIDs( TALEFamily[] families, TALEFamilyBuilder builder ) {
			for(int i=0,k=1;i<families.length;i++){
				if(families[i].getFamilyId() == null){
					families[i].setFamilyId( "New Class "+(k) );
				}
			}			
		}
		
		
		
	}
	
	public interface FamilyIdGenerator{
		
		public void setFamilyIDs(TALEFamily[] families, TALEFamilyBuilder builder);
		
	}
	
	public static class TALEFamily implements PlotGenerator, Comparable<TALEFamily> {
		
		private String id;
		private StringAlignment[][] alignments;
		private ClusterTree<TALE> tree;
		
		
		
		private TALEFamily(String familyId, ClusterTree<TALE> tree, TALEFamilyBuilder builder){
			this.tree = tree;
			this.id = familyId;
		/*	this.costs = costs;
			this.linkage = linkage;
			this.at = at;
			this.cut = cut;
			this.extraGapOpening = extraGapOpening;
			this.extraGapExtension = extraGapExtension;*/
			
			TALE[] members = tree.getClusterElements();
			alignments = new StringAlignment[members.length][members.length];
			for(int j=0;j<members.length;j++){
				for(int k=0;k<members.length;k++){
					if(j != k){
						StringAlignment sa = TALEAligner.align( members[j], members[k], builder.costs, builder.at, builder.extraGapOpening, builder.extraGapExtension );
						alignments[j][k] = sa;
					}
				}
			}
			
		}
		
		public TALEFamily(StringBuffer xml) throws NonParsableException{
			xml = XMLParser.extractForTag(xml,"TALEClass");
			alignments = (StringAlignment[][])XMLParser.extractObjectForTags( xml, "alignments" );
			/*at = (AlignmentType)XMLParser.extractObjectForTags( xml, "at" );
			costs = (Costs)XMLParser.extractObjectForTags( xml, "costs" );
			cut = (Double)XMLParser.extractObjectForTags( xml, "cut" );
			extraGapExtension = (double)XMLParser.extractObjectForTags( xml, "extraGapExtension" );
			extraGapOpening = (double)XMLParser.extractObjectForTags( xml, "extraGapOpening" );*/
			id = (String)XMLParser.extractObjectForTags( xml, "id" );
			//linkage = (Linkage)XMLParser.extractObjectForTags( xml, "linkage" );
			tree = (ClusterTree<TALE>)XMLParser.extractObjectForTags( xml, "tree" );
			
		}
		
		
		public int compareTo(TALEFamily tf2){
			return id.compareTo( tf2.getFamilyId() );
		}
		
		public StringBuffer toXML(){
			StringBuffer xml = new StringBuffer();
			XMLParser.appendObjectWithTags( xml, alignments, "alignments" );
/*			XMLParser.appendObjectWithTags( xml, at, "at" );
			XMLParser.appendObjectWithTags( xml, costs, "costs" );
			XMLParser.appendObjectWithTags( xml, cut, "cut" );
			XMLParser.appendObjectWithTags( xml, extraGapExtension, "extraGapExtension" );
			XMLParser.appendObjectWithTags( xml, extraGapOpening, "extraGapOpening" );*/
			XMLParser.appendObjectWithTags( xml, id, "id" );
			//XMLParser.appendObjectWithTags( xml, linkage, "linkage" );
			XMLParser.appendObjectWithTags( xml, tree, "tree" );
			XMLParser.addTags( xml, "TALEClass" );
			return xml;
		}
		
		
		public String getSpecificityConsensus(TALgetterDiffSM model) throws IllegalArgumentException, WrongAlphabetException{
			double[][] spec = getSpecificityProfile( model );
			AlphabetContainer dna = DNAAlphabetContainer.SINGLETON;
			StringBuffer sb = new StringBuffer();
			sb.append( dna.getSymbol( 0, ToolBox.getMaxIndex( spec[0] ) ) );
			for(int i=1;i<spec.length;i++){
				sb.append( dna.getSymbol( 0, ToolBox.getMaxIndex( spec[i] ) ) );
				if(i<spec.length-1){
					sb.append( "  " );
				}
			}
			return sb.toString();
		}
		
		public double[][] getSpecificityProfile(TALgetterDiffSM model) throws IllegalArgumentException, WrongAlphabetException{
			
			Pair<TALE[],String[]> al = this.getInducedMultipleAlignment();
			
			DiscreteAlphabet modelCon = (DiscreteAlphabet)model.getRVDAlphabet().getAlphabetAt( 0 );
			
			String[] als = al.getSecondElement();
			double[][] specs = null;
			for(int i=0;i<als.length;i++){
				String[] parts = als[i].trim().split( " " );
				if(specs == null){
					specs = new double[parts.length+1][4];
				}
				double[][] temp = new double[parts.length+1][4];
				for(int j=0;j<parts.length;j++){
					
					if(parts[j].equals( "--" ) || !modelCon.isSymbol( parts[j] )){
						temp[j+1] = new double[4];
						Arrays.fill( temp[j+1], 0.25 );
					}else{
						Sequence seq = Sequence.create( model.getRVDAlphabet(), parts[j] );
						Pair<double[][],double[]> pair = model.getSpecificitiesAndImportances( seq );
						temp[j+1] = pair.getFirstElement()[1];
						for(int k=0;k<temp[j+1].length;k++){
							temp[j+1][k] = temp[j+1][k]*pair.getSecondElement()[0] + 0.25*(1.0-pair.getSecondElement()[0]);
						}
						if(j==0 || parts[j-1].equals( "--" )){
							temp[j] = pair.getFirstElement()[0];
						}
					}
				}
				for(int j=0;j<temp.length;j++){
					for(int k=0;k<temp[j].length;k++){
						specs[j][k] += temp[j][k]/(double)als.length;
					}
				}
			}
			return specs;
		}
		
		
		public String getFamilyId(){
			return id;
		}
		
		public void setFamilyId(String newId){
			this.id = newId;
		}
		
		public StringAlignment getAlignmentForIDs(String id1, String id2){
			TALE[] members = getFamilyMembers();
			for(int i=0;i<members.length;i++){
				if(id1.equals( members[i].getId() )){
					
					for(int j=0;j<members.length;j++){
						if(id2.equals( members[j].getId() )){
							return alignments[i][j];
						}
					}
					
				}
			}
			//System.out.println("did not find alignment for "+id1+"; "+id2);
			//System.out.println("members: "+Arrays.toString( members ));
			return null;
		}
		
		private FamilyDistance getDist(TALEFamilyBuilder builder){
			FamilyDistance dist = null;
			if(builder.linkage == Linkage.SINGLE){
				dist = FamilyDistance.MIN;
			}else if(builder.linkage == Linkage.COMPLETE){
				dist = FamilyDistance.MAX;
			}else if(builder.linkage == Linkage.AVERAGE){
				dist = FamilyDistance.MEAN;
			}
			return dist;
		}
		
		public double getDistance(TALE tale, FamilyDistance dist, TALEFamilyBuilder builder){
			if(dist == null){
				dist = getDist(builder);
			}
			TALE[] members = getFamilyMembers();
			double[] ds = new double[members.length];
			for(int i=0;i<members.length;i++){
				ds[i] = TALEAligner.align( tale, members[i], builder.costs, builder.at, builder.extraGapOpening, builder.extraGapExtension ).getCost();
			}
			if(dist == FamilyDistance.MAX){
				return ToolBox.max( ds );
			}else if(dist == FamilyDistance.MIN){
				return ToolBox.min( ds );
			}else if(dist == FamilyDistance.MEAN){
				return ToolBox.mean( 0, ds.length, ds );
			}else{
				return Double.NaN;
			}
		}
		
		
		public double getSignificance(TALE tale, AlignmentPValues pv, FamilyDistance dist, TALEFamilyBuilder builder){
			if(dist == null){
				dist = getDist(builder);
			}
			if(pv == null){
				if(builder.costs instanceof AffineCosts){
					Costs c2 = ((AffineCosts)builder.costs).getInternalCosts();
					if(c2 instanceof RVDCosts){
						RVDCosts c3 = (RVDCosts) c2;
						pv = new AlignmentPValues( builder.getAllTALEs(), c3 );
					}
				}
			}
			TALE[] members = getFamilyMembers();
			double[] ds = new double[members.length];
			for(int i=0;i<members.length;i++){
				double temp = TALEAligner.align( tale, members[i], builder.costs, builder.at, builder.extraGapOpening, builder.extraGapExtension ).getCost();
				ds[i] = pv.getLog10PValue( tale, members[i], temp, builder.extraGapOpening, builder.extraGapExtension );
			}
			if(dist == FamilyDistance.MAX){
				return ToolBox.max( ds );
			}else if(dist == FamilyDistance.MIN){
				return ToolBox.min( ds );
			}else if(dist == FamilyDistance.MEAN){
				double di = 0;
				for(int i=0;i<ds.length;i++){
					di += Math.log1p( -Math.exp(ds[i]*Math.log( 10 )) );
				}
				
				return Math.log1p( -Math.exp(di) )/Math.log( 10 );
			}else{
				return Double.NaN;
			}
		}
		
		
		public double getFamilySignificance(AlignmentPValues pv, TALEFamilyBuilder builder){
			
			if(pv == null){
				if(builder.costs instanceof AffineCosts){
					Costs c2 = ((AffineCosts)builder.costs).getInternalCosts();
					if(c2 instanceof RVDCosts){
						RVDCosts c3 = (RVDCosts) c2;
						pv = new AlignmentPValues( builder.getAllTALEs(), c3 );
					}
				}
			}
			
			if(pv == null){
				return Double.NaN;
			}
			
			double famsig = 0.0;
			
			TALE[] members = getFamilyMembers();
			
			for(int i=1;i<alignments.length;i++){
				for(int j=0;j<i;j++){
					double p = pv.getLog10PValue( members[i], members[j], alignments[i][j].getCost(), builder.extraGapOpening, builder.extraGapExtension );
					p = Math.log1p( -Math.exp(p*Math.log( 10 )) );
					famsig += p;
				}
			}
			return Math.log1p( -Math.exp( famsig ) )/Math.log( 10 );
		}
		
		
		public TALE[] getFamilyMembers(){
			return tree.getClusterElements();
		}
		
		public int getFamilySize(){
			return tree.getNumberOfElements();
		}
		
		public String toString(TALgetterDiffSM model, TALEFamilyBuilder builder){
			StringBuffer sb = new StringBuffer();
			sb.append("Class "+id+" for ("+builder.costs.getClass().getSimpleName()+", "+builder.cut+", "+builder.at+")\n");
			sb.append( "distance: "+format.format( tree.getDistance())+"\n" );
			
			if(builder.costs instanceof AffineCosts){
				Costs c2 = ((AffineCosts)builder.costs).getInternalCosts();
				if(c2 instanceof RVDCosts){
					RVDCosts c3 = (RVDCosts) c2;
					double p = this.getFamilySignificance( new AlignmentPValues( builder.getAllTALEs(), c3 ), builder );
					sb.append( "significance: p="+formatE.format( Math.pow( 10, p ) )+"\n" );
				}
			}
			sb.append( "\n"+this.inducedMultipleAlignmentToString()+"\n" );
			if(model != null){
				try{
					String cons = getSpecificityConsensus( model );
					sb.append( "Most likely common binding sequence:\n" );
					sb.append( cons );
					sb.append( "\n" );
				}catch(WrongAlphabetException e){ }
			}
			sb.append( "\n\nClass tree:\n" );
			sb.append(tree.toNewick()+"\n\n");
			sb.append( "Alignment scores:\n" );
			TALE[] members = tree.getClusterElements();
			for(int i=1;i<members.length;i++){
				for(int j=0;j<i;j++){
					sb.append( members[i].getId()+" vs. "+members[j].getId()+": "+format.format( alignments[i][j].getCost() )+"\n" );
				}
			}
			
			return sb.toString();
		}

		public ClusterTree<TALE> getTree() {
			return tree;
		}
		
		public void generatePlot(GraphicsAdaptor adaptor) throws IOException {
			Graphics2D dummy = adaptor.getGraphics( 10, 10 );
			
			TALEFamilyTreePlotter plotter = new TALEFamilyTreePlotter( 30 );
			
			int[] dim = plotter.getDimension( dummy, this );
			
			
			Graphics2D graphics = adaptor.getGraphics( dim[0], dim[1] );
			
			graphics.setColor( Color.white );
			
			graphics.fillRect( 0, 0, dim[0], dim[1] );
			
			graphics.setColor( Color.black );
			
			graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			
			plotter.plot( graphics, this );
		}
		
		public void plotFamilyToFile(String nameBase, GraphicsAdaptor adaptor) throws IOException{
			
			generatePlot( adaptor );
		
			adaptor.generateOutput( nameBase+"."+adaptor.getGraphicsExtension() );
			
		}

		
		private TALEFamily removeTALE( double[][] dmat,  LinkedList<ClusterTree<TALE>> tales, int indexOff, TALEFamilyBuilder builder){
			ClusterTree<TALE>[] leaves = this.tree.getLeaves();
			LinkedList<ClusterTree<TALE>> list = new LinkedList<ClusterTree<TALE>>();
			Collections.addAll( list, leaves );
			System.out.println("before "+list);
			System.out.println("to remove "+tales);
			
			list.removeAll(tales);
			System.out.println("after "+list);
			if(list.size() == 0){
				return null;
			}else{
				Hclust<TALE> hclust = new Hclust<TALE>( null, builder.linkage );
				
				ClusterTree<TALE> newTree = hclust.cluster( indexOff, dmat, list.toArray( new ClusterTree[0] ) );//TODO efficiency
				
			//	System.out.println(newTree);
				
				//ClusterTree<TALE> newTree = cluster(newTales,linkage,costs,at, extraGapOpening, extraGapExtension).getSecondElement();
				
				newTree.leafOrder(dmat);
				
				TALEFamilyBuilder.renameTALEs(newTree, this.getFamilyId());
				
				TALEFamily fam = new TALEFamily( this.getFamilyId(), newTree, builder );
				
				return fam;
			}
		}
		
		private TALEFamily addTALE( double[][] dmat, TALE tale, int newIndex, int indexOff, TALEFamilyBuilder builder ) {
			
			ClusterTree<TALE>[] leaves = this.tree.getLeaves();
			LinkedList<ClusterTree<TALE>> list = new LinkedList<ClusterTree<TALE>>();
			Collections.addAll( list, leaves );
			list.add( new ClusterTree<TALE>( tale, newIndex ) );
		//	System.out.println("Adding to family "+id);
			//for(int i=0;i<list.size();i++){
			//	System.out.println(i+" "+list.get(i).getClusterElements()[0]+" "+list.get( i ).getOriginalIndex());
				//for(int j=0;j<i;j++){
				//	System.out.println("    "+j+" "+list.get( j ).getOriginalIndex());
				//	System.out.println(dmat[list.get( i ).getOriginalIndex()][list.get( j ).getOriginalIndex()]);
				//}
			//}
			
			
			Hclust<TALE> hclust = new Hclust<TALE>( null, builder.linkage );
			
			ClusterTree<TALE> newTree = hclust.cluster( indexOff, dmat, list.toArray( new ClusterTree[0] ) );//TODO efficiency
			
		//	System.out.println(newTree);
			
			//ClusterTree<TALE> newTree = cluster(newTales,linkage,costs,at, extraGapOpening, extraGapExtension).getSecondElement();
			
			newTree.leafOrder(dmat);
			TALEFamily fam = new TALEFamily( this.getFamilyId(), newTree, builder );
			
			return fam;
			
		}
		
		
		public String inducedMultipleAlignmentToString(){
			
			StringBuffer sb = new StringBuffer();
			
			Pair<TALE[],String[]> al = getInducedMultipleAlignment();
			
			String[] als = al.getSecondElement();
			TALE[] tales = al.getFirstElement();
			for(int i=0;i<als.length;i++){
				sb.append(als[i]+"\t"+tales[i].getId()+"\n");
			}
			
			return sb.toString();
		}
		
		public Pair<TALE[],String[]> getInducedMultipleAlignment(){
			return getInducedMultipleAlignment( getTree(), this );
		}
		
		
		private Pair<TALE[],String[]> getInducedMultipleAlignment2(ClusterTree<TALE> tree, TALEFamily family){
			TALE[] members = tree.getClusterElements();
			if(members.length == 1){
				Sequence rvds = members[0].getRvdSequence();
				return new Pair<TALE[],String[]>(new TALE[]{members[0]},new String[]{" "+rvds.toString() });
			}else if(members.length == 2){
				//System.out.println("searching alignment for IDs "+members[0].getId()+","+members[1].getId());
				StringAlignment al = family.getAlignmentForIDs( members[0].getId(), members[1].getId() );
				return new Pair<TALE[],String[]>( new TALE[]{members[0],members[1]},new String[]{al.getAlignedString( 0 ),al.getAlignedString( 1 )} );
			}else{
				
				ClusterTree<TALE>[] subs = tree.getSubTrees();
				
				TALE top = subs[0].getClusterElements()[subs[0].getNumberOfElements()-1];
				TALE bot = subs[1].getClusterElements()[0];
				
				Pair<TALE[],String[]> topPair = getInducedMultipleAlignment( subs[0], family );
				Pair<TALE[],String[]> botPair = getInducedMultipleAlignment( subs[1], family );
				
				String[] topTemp = topPair.getSecondElement();
				StringBuffer[] topAls = new StringBuffer[topTemp.length];
				for(int i=0;i<topAls.length;i++){
					topAls[i] = new StringBuffer(topTemp[i].replaceAll(" ", ""));
				}
				String[] botTemp = botPair.getSecondElement();
				StringBuffer[] botAls = new StringBuffer[botTemp.length];
				for(int i=0;i<botAls.length;i++){
					botAls[i] = new StringBuffer(botTemp[i].replaceAll(" ", ""));
				}
				
				StringAlignment al = family.getAlignmentForIDs( top.getId(), bot.getId() );
				
				String al2 = al.getAlignedString( 0 ).replaceAll( " ", "" );
				String al3 = al.getAlignedString( 1 ).replaceAll( " ", "" );
				
				int frontTopNew = 0;
				while( al2.charAt(frontTopNew) == '-' ){ frontTopNew++; }
				int frontTopOld = 0;
				while(topAls[topAls.length-1].charAt(frontTopOld) == '-'){ frontTopOld++; }
				int frontBotNew = 0;
				while( al3.charAt(frontBotNew) == '-' ){ frontBotNew++; }
				int frontBotOld = 0;
				while( botAls[0].charAt(frontBotOld) == '-' ){ frontBotOld++; }
				
				int frontTop = Math.max(frontTopNew, frontTopOld);
				int frontBot = Math.max(frontBotNew, frontBotOld);
				

				
				for(int i=frontTopOld;i<frontTop;i++){
					for(int j=0;j<topAls.length;j++){
						topAls[j].insert(0, "-");
					}
				}
				
				for(int i=frontBotOld;i<frontBot;i++){
					for(int j=0;j<botAls.length;j++){
						botAls[j].insert(0, "-");
					}
				}
				
				
				
				int maxLen = 0;
				
				for(int i=0;i<topAls.length;i++){
					if(topAls[i].length() > maxLen){
						maxLen = topAls[i].length(); 
					}
				}
				
				for(int i=0;i<botAls.length;i++){
					if(botAls[i].length() > maxLen){
						maxLen = botAls[i].length();
					}
				}
				
				for(int i=0;i<topAls.length;i++){
					for(int j=topAls[i].length();j<maxLen;j++){
						topAls[i].append("-");
					}
				}
				
				for(int i=0;i<botAls.length;i++){
					for(int j=botAls[i].length();j<maxLen;j++){
						botAls[i].append("-");
					}
				}
				
				
				
				LinkedList<TALE> resTales = new LinkedList<TALE>();
				LinkedList<String> resAls = new LinkedList<String>();
				
				Collections.addAll(resTales, topPair.getFirstElement());
				Collections.addAll(resTales, botPair.getFirstElement());
				
				for(int i=0;i<topAls.length;i++){
					for(int k=topAls[i].length()-2;k>=0;k-=2){
						topAls[i].insert(k, " ");
					}
					resAls.add(topAls[i].toString());
				}
				for(int i=0;i<botAls.length;i++){
					for(int k=botAls[i].length()-2;k>=0;k-=2){
						botAls[i].insert(k, " ");
					}
					resAls.add(botAls[i].toString());
				}
				
				return new Pair<TALE[],String[]>(resTales.toArray( new TALE[0] ), resAls.toArray( new String[0] ));
				
			}
			
			
		}
		
		
		private Pair<TALE[],String[]> getInducedMultipleAlignment(ClusterTree<TALE> tree, TALEFamily family){
			TALE[] members = tree.getClusterElements();
			if(members.length == 1){
				Sequence rvds = members[0].getRvdSequence();
				return new Pair<TALE[],String[]>(new TALE[]{members[0]},new String[]{" "+rvds.toString() });
			}else if(members.length == 2){
				//System.out.println("searching alignment for IDs "+members[0].getId()+","+members[1].getId());
				StringAlignment al = family.getAlignmentForIDs( members[0].getId(), members[1].getId() );
				return new Pair<TALE[],String[]>( new TALE[]{members[0],members[1]},new String[]{al.getAlignedString( 0 ),al.getAlignedString( 1 )} );
			}else{
				ClusterTree<TALE>[] subs = tree.getSubTrees();
				TALE[][] tales = new TALE[subs.length][];
				StringBuffer[][] als = new StringBuffer[subs.length][];
				for(int i=0;i<subs.length;i++){
					Pair<TALE[],String[]> pair = getInducedMultipleAlignment( subs[i], family );
					tales[i] = pair.getFirstElement();
					String[] temp = pair.getSecondElement();
					als[i] = new StringBuffer[temp.length];
					for(int j=0;j<als[i].length;j++){
						als[i][j] = new StringBuffer(temp[j].replaceAll( " ", "" ));
					}
				}
				
				TALE prev = tales[0][tales[0].length-1];
				String prevAl = als[0][als[0].length-1].toString();
				
				
				for(int i=1;i<subs.length;i++){
					TALE curr = tales[i][0];
					String currAl = als[i][0].toString();
					
					StringAlignment al = family.getAlignmentForIDs( prev.getId(), curr.getId() );
					
					
					
					String al2 = al.getAlignedString( 0 ).replaceAll( " ", "" );
					String al3 = al.getAlignedString( 1 ).replaceAll( " ", "" );
					
					
					
					//Copy gaps from original top alignment to new alignment and vice versa
					
					int idx1 = prevAl.length()-1, idx2 = al2.length()-1;
					int[] addPrev = new int[prevAl.length()+1];
					int[] addCurr = new int[al2.length()+1];
					while( idx1 >=0 || idx2 >=0 ){
						if(idx1 >= 0 && idx2 >= 0 && prevAl.charAt( idx1 ) == al2.charAt( idx2 )){//indentical
							idx1--;
							idx2--;
						}else if(idx1 >= 0 && prevAl.charAt( idx1 ) == '-'){//we had a gap that we don't have now
							addCurr[idx2+1]++;
							idx1--;
						}else if(idx2 >= 0 && al2.charAt( idx2 ) == '-'){
							addPrev[idx1+1]++;
							idx2--;
						}
						
					}
					
					//insert gaps into top original alignment
					for(int j=0;j<i;j++){
						for(int k=0;k<als[j].length;k++){
							for(int l=addPrev.length-1;l>=0;l--){
								for(int n=0;n<addPrev[l];n++){
									als[j][k].insert( l, '-' );
								}
							}
						}
					}
					
					
					//insert gaps into second aligned string as well
					StringBuffer temp = new StringBuffer(al3);
					for(int l=addCurr.length-1;l>=0;l--){
						for(int n=0;n<addCurr[l];n++){
							temp.insert( l, '-' );
						}
					}
					al3 = temp.toString();
					
					
									
					//copy (new or old) gaps from second aligned string to original bottom alignment
					//and gaps in original bottom alignment to second aligned string *and* original top alignment
					idx1 = currAl.length()-1; idx2 = al3.length()-1;
					addCurr = new int[currAl.length()+1];
					addPrev = new int[al3.length()+1];
					while( idx1 >= 0 || idx2 >= 0 ){
						if(idx1 >= 0 && idx2 >= 0 && currAl.charAt(idx1) == al3.charAt( idx2 )){
							idx1--;
							idx2--;
						}else if(idx1>=0 && currAl.charAt( idx1 ) == '-'){
							addPrev[idx2+1]++;
							idx1--;
						}else if(idx2>=0 && al3.charAt( idx2 ) == '-'){
							addCurr[idx1+1]++;
							idx2--;
						}
					}
				
					//insert gaps into original top alignment
					for(int j=0;j<i;j++){
						for(int k=0;k<als[j].length;k++){
							for(int l=addPrev.length-1;l>=0;l--){
								for(int n=0;n<addPrev[l];n++){
									als[j][k].insert( l, '-' );
								}
							}
						}
						
					}
					
					//insert gaps into original bottom alignment
					for(int k=0;k<als[i].length;k++){
						for(int l=addCurr.length-1;l>=0;l--){
							for(int n=0;n<addCurr[l];n++){
								als[i][k].insert( l, '-' );
							}
						}
					}
					
					prev = tales[i][tales[i].length-1];
					prevAl = als[i][als[i].length-1].toString();
				}
				
				
				/**
				 * from here
				 */
				LinkedList<StringBuffer> tempAls = new LinkedList<StringBuffer>();
				for(int i=0;i<tales.length;i++){
					for(int j=0;j<tales[i].length;j++){
						tempAls.add( als[i][j] );
					}
				}
				
				
				boolean[][] onlyGaps = new boolean[tempAls.size()][tempAls.get(0).length()];
				boolean[][] onlyGapsEnd = new boolean[tempAls.size()][tempAls.get(0).length()];
				for(int i=0;i<onlyGaps.length;i++){
					int j=0;
					while(tempAls.get(i).charAt(j)=='-'){ onlyGaps[i][j] = true;j++; }
					j=tempAls.get(i).length()-1;
					while(tempAls.get(i).charAt(j)=='-'){ onlyGapsEnd[i][j] = true;j--; }
				}
				for(int i=0;i<onlyGaps[0].length;i++){
					int j=0;
					while(j<tempAls.size() && tempAls.get(j).charAt(i)=='-' && !onlyGaps[j][i]){j++;}
					if(j>0 && j<tempAls.size() && tempAls.get(j).charAt(i)=='-'){
						for(int k=0;k<j;k++){
							tempAls.get(k).delete(i, i+1);
							tempAls.get(k).insert(0, '-');
						}
					}
					j=tempAls.size()-1;
					while(j>=0 && tempAls.get(j).charAt(i)=='-' && !onlyGaps[j][i]){j--;}
					if(j>0 && j<tempAls.size() && tempAls.get(j).charAt(i)=='-'){
						for(int k=tempAls.size()-1;k>j;k--){
							tempAls.get(k).delete(i, i+1);
							tempAls.get(k).insert(0, '-');
						}
					}
				}
				onlyGaps = onlyGapsEnd;
				for(int i=onlyGaps[0].length-1;i>=0;i--){
					int j=0;
					while(j<tempAls.size() && tempAls.get(j).charAt(i)=='-' && !onlyGaps[j][i]){j++;}
					if(j>0 && j<tempAls.size() && tempAls.get(j).charAt(i)=='-'){
						for(int k=0;k<j;k++){
							tempAls.get(k).delete(i, i+1);
							tempAls.get(k).append('-');
						}
					}
					j=tempAls.size()-1;
					while(j>=0 && tempAls.get(j).charAt(i)=='-' && !onlyGaps[j][i]){j--;}
					if(j>0 && j<tempAls.size() && tempAls.get(j).charAt(i)=='-'){
						for(int k=tempAls.size()-1;k>j;k--){
							tempAls.get(k).delete(i, i+1);
							tempAls.get(k).append('-');
						}
					}
				}
				
					
					
				
				LinkedList<TALE> resTales = new LinkedList<TALE>();
				LinkedList<String> resAls = new LinkedList<String>();
				
				for(int i=0;i<tales.length;i++){
					for(int j=0;j<tales[i].length;j++){
						resTales.add( tales[i][j] );
					}
				}
				for(int i=0;i<tempAls.size();i++){
					for(int k=tempAls.get(i).length()-2;k>=0;k-=2){
						tempAls.get(i).insert( k, ' ' );
					}
					resAls.add(tempAls.get(i).toString());
				}
				
				/**
				 * to here
				 */
				
				/*for(int i=0;i<tales.length;i++){
					for(int j=0;j<tales[i].length;j++){
						resTales.add( tales[i][j] );
						for(int k=als[i][j].length()-2;k>=0;k-=2){
							als[i][j].insert( k, ' ' );
						}
						resAls.add( als[i][j].toString() );
					}
				}*/
				
				return new Pair<TALE[],String[]>(resTales.toArray( new TALE[0] ), resAls.toArray( new String[0] ));
			}
			
		}
		
		
	} 
	
	private TALEFamily[] families;
	private double[][] dmat;
	//private ClusterTree<TALE> familyTree;
	
	private Costs costs;
	private AlignmentType at;
	

	private Linkage linkage;
	private double extraGapOpening;
	private double extraGapExtension;
	private double cut;
	private double pval;
	private String[] reservedNames;
	
	public TALEFamilyBuilder(TALE[] tales) throws IllegalArgumentException, IOException, WrongAlphabetException{
		this(tales, new AffineCosts(5.0, new RVDCosts( 1.0, 0.2, 0.8, 0.0 ) ), Linkage.AVERAGE, AlignmentType.SEMI_GLOBAL, 1.0,0.1, 5.0, 0.01 );
	}
	
	public TALEFamilyBuilder(TALE[] tales, Costs costs, Linkage linkage, AlignmentType at, double extraGapOpening, double extraGapExtension, double cut, double pval) {
		this.pval = pval;
		this.at = at;
		this.costs = costs;
		this.cut = cut;
		this.extraGapExtension = extraGapExtension;
		this.extraGapOpening = extraGapOpening;
		this.linkage = linkage;
		
		Pair<double[][],ClusterTree<TALE>> pair = cluster(tales,linkage,costs,at, extraGapOpening, extraGapExtension);
		
		ClusterTree<TALE> tree = pair.getSecondElement();
		
		//this.familyTree = tree;
		
		dmat = pair.getFirstElement();
		
		
		ClusterTree<TALE>[] subtrees = Hclust.cutTree( cut, tree );	
		
		
		this.families = new TALEFamily[subtrees.length];
		for(int i=0;i<subtrees.length;i++){
			
			subtrees[i].leafOrder( dmat );
			this.families[i] = new TALEFamily( (i+1)+"", subtrees[i], this );
		}
		
	}
	
	public TALEFamilyBuilder(StringBuffer xml) throws NonParsableException {
		at = (AlignmentType)XMLParser.extractObjectForTags( xml, "at" );
		costs = (Costs)XMLParser.extractObjectForTags( xml, "costs" );
		cut = (Double)XMLParser.extractObjectForTags( xml, "cut" );
		StringBuffer sb = XMLParser.extractForTag(xml, "dmatStore" );
		if(sb != null){
			dmat = parseDmat(sb);
		}else{
			dmat = (double[][])XMLParser.extractObjectForTags( xml, "dmat" );
		}
		extraGapOpening = (Double)XMLParser.extractObjectForTags( xml, "extraGapOpening" );
		extraGapExtension = (Double)XMLParser.extractObjectForTags( xml, "extraGapExtension" );
		StringBuffer tempXML = new StringBuffer(xml);
		try{
			families = (TALEFamily[]) XMLParser.extractObjectForTags(xml, "families");
		}catch(Exception ex){
			xml = tempXML;
			String[] temp = (String[])XMLParser.extractObjectForTags( xml, "families" );
			families = new TALEFamily[temp.length];
			for(int i=0;i<families.length;i++){
				families[i] = new TALEFamily( new StringBuffer( temp[i] ) );
			}
		}
		linkage = (Linkage)XMLParser.extractObjectForTags( xml, "linkage" );
		reservedNames = (String[])XMLParser.extractObjectForTags( xml, "reservedNames" );
		pval = (Double)XMLParser.extractObjectForTags( xml, "pval" );
	}
	
	public AlignmentType getAlignmentType() {
		return at;
	}
	
	public double getExtraGapOpening() {
		return extraGapOpening;
	}

	public double getExtraGapExtension() {
		return extraGapExtension;
	}

	public StringBuffer toXML(){
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, at, "at" );
		XMLParser.appendObjectWithTags( xml, costs, "costs" );
		XMLParser.appendObjectWithTags( xml, cut, "cut" );
		//XMLParser.appendObjectWithTags( xml, dmat, "dmat" );
		StringBuffer dmat2 = storeDmat(dmat);
		XMLParser.addTags(dmat2, "dmatStore");
		xml.append(dmat2);
		XMLParser.appendObjectWithTags( xml, extraGapOpening, "extraGapOpening" );
		XMLParser.appendObjectWithTags( xml, extraGapExtension, "extraGapExtension" );
		/*String[] temp = new String[families.length];
		for(int i=0;i<temp.length;i++){
			temp[i] = families[i].toXML().toString();
		}*/
		XMLParser.appendObjectWithTags( xml, families, "families" );
		XMLParser.appendObjectWithTags( xml, linkage, "linkage" );
		XMLParser.appendObjectWithTags( xml, reservedNames, "reservedNames" );
		XMLParser.appendObjectWithTags( xml, pval, "pval" );
		return xml;
	}
	
	private static StringBuffer storeDmat(double[][] dmat){
		StringBuffer sb = new StringBuffer();
		DecimalFormat df = new DecimalFormat("#.######");
		for(int i=0;i<dmat.length;i++){
			for(int j=0;j<dmat[i].length;j++){
				if(j>0){
					sb.append(";");
				}
				sb.append(df.format(dmat[i][j]));
			}
			sb.append("$");
		}
		return sb;
	}
	
	
	private static double[][] parseDmat(StringBuffer sb){
		int off = 0;
		int n = 0;
		while( ( off = sb.indexOf("$",off)+1 ) > 0){
			n++;
		}
		double[][] dmat = new double[n][];
		
		off = 0;
		int off2 = 0;
		n=0;
		while( (off2 = sb.indexOf("$", off)) >= 0 ){
			
			String[] parts = sb.substring(off, off2).split(";");
			dmat[n] = new double[parts.length];
			for(int j=0;j<parts.length;j++){
				parts[j] = parts[j].replaceAll(",", ".");
				dmat[n][j] = Double.parseDouble(parts[j]);
			}
			n++;
			off = off2+1;
		}
		return dmat;
	}
	
	
	
/*	public ClusterTree<TALE> getFamilyTree(){
		return familyTree;
	}
	
	public double[][] getDistanceMatrix() throws CloneNotSupportedException{
		int numEls = familyTree.getNumberOfElements();
		double[][] dmat2 = new double[numEls][numEls];
		ClusterTree<TALE>[] leaves = familyTree.getLeaves();
		for(int i=0;i<leaves.length;i++){
			for(int j=0;j<leaves.length;j++){
				dmat2[i][j] = dmat[ leaves[i].getOriginalIndex() ][ leaves[j].getOriginalIndex() ];
			}
		}
		return dmat2;
	}*/
	
	
	/*public TALEFamily[] rebuild(Linkage linkage, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		
		LinkedList<ClusterTree<TALE>> leaves = new LinkedList<ClusterTree<TALE>>();
		int minIndex = 0;
		for(int i=0;i<families.length;i++){
			ClusterTree<TALE> famTree = families[i].getTree();
			int index = famTree.getMinimumOriginalIndex();
			if(index < minIndex){
				minIndex = index;
			}
			ClusterTree<TALE>[] myLeaves = famTree.getLeaves();
			for(int j=0;j<myLeaves.length;j++){
				leaves.add( myLeaves[j] );
			}
		}
		
		TALE[] tales = new TALE[leaves.size()];
		for(int i=0;i<leaves.size();i++){
			ClusterTree<TALE> leaf = leaves.get( i );
			tales[leaf.getOriginalIndex()] = leaf.getClusterElements()[0];
		}
		
		double[][] dmat = computeDistMatrix( tales, costs, at, extraGapOpening, extraGapExtension );
		Hclust<TALE> clust = new Hclust<>( null, linkage );
		
		TALEFamily[] fams = new TALEFamily[families.length];
		for(int i=0;i<families.length;i++){
			ClusterTree<Integer> intTree = families[i].getTree().getIndexTree();
			ClusterTree<Integer>[] myLeaves = intTree.getLeaves();
			LinkedList<ClusterTree<Integer>> leafList = new LinkedList<ClusterTree<Integer>>();
			for(int j=0;j<myLeaves.length;j++){
				leafList.add( myLeaves[j] );
			}
			
			ClusterTree<Integer> reclustered = clust.cluster( dmat, leafList, -minIndex );
			if(!reclustered.isLeaf()){
				reclustered = new ClusterTree<Integer>(reclustered.getDistance(),intTree.getOriginalIndex(),reclustered.getSubTrees());
				reclustered.leafOrder( this.dmat );//TODO???
				int index = reclustered.getMinimumOriginalIndex();
				if(index < minIndex){
					minIndex = index;
				}
			}
			ClusterTree<TALE> re2 = clust.createTree( reclustered, tales );
			fams[i] = new TALEFamily( families[i].getFamilyId(), re2, costs, at, re2.getDistance(), linkage, extraGapOpening, extraGapExtension );
		}
		
		/*clust = new Hclust<>( null, Linkage.SINGLE );
		
		for(int i=0;i<fams.length;i++){
			double[] ds = new double[fams.length];
			ds[i] = Double.POSITIVE_INFINITY;
			for(int j=0;j<fams.length;j++){
				if(i != j){
					ds[j] = clust.getDistance( dmat, fams[i].getTree().getIndexTree(), fams[j].getTree().getIndexTree() );
				}
			}
			System.out.println("fam "+(i+1)+"\t"+fams[i].getTree().getDistance()+"\t"+ToolBox.min( ds ));
			
			
		}
		
		
		
		return fams;
	}
	*/
	
	public double getCut(){
		return cut;
	}
	
	public void setToOld(){
		for(int i=0;i<families.length;i++){
			TALE[] mems = families[i].getFamilyMembers();
			for(int j=0;j<mems.length;j++){
				mems[j].setIsNew( false );
			}
		}
	}
	
	public TALE[] getAllTALEs(){
		LinkedList<TALE> all = new LinkedList<TALE>();
		for(int i=0;i<families.length;i++){
			TALE[] temp = families[i].getFamilyMembers();
			for(int j=0;j<temp.length;j++){
				all.add( temp[j] );
			}
		}
		return all.toArray( new TALE[0] );
	}
	
	public ClusterTree<TALEFamily> clusterFamilies(){
		
		LinkedList<ClusterTree<Integer>> intTrees = new LinkedList<ClusterTree<Integer>>();
		LinkedList<ClusterTree<TALE>> leaves = new LinkedList<ClusterTree<TALE>>();
		IntList familyRootOriginalIndexes = new IntList();
		int minIndex = 0;
		for(int i=0;i<families.length;i++){
			ClusterTree<TALE> famTree = families[i].getTree();
			int index = famTree.getMinimumOriginalIndex();
			if(index < minIndex){
				minIndex = index;
			}
			familyRootOriginalIndexes.add( famTree.getOriginalIndex() );
			//System.out.println(i+" "+families[i].getFamilyId()+" "+famTree.getOriginalIndex());
			//TALE[] tales = famTree.getClusterElements();
			ClusterTree<Integer> indexTree = famTree.getIndexTree( );
			ClusterTree<TALE>[] myLeaves = famTree.getLeaves();
			for(int j=0;j<myLeaves.length;j++){
				leaves.add( myLeaves[j] );
			}
			intTrees.add( indexTree );
		}
		
		TALE[] tales = new TALE[leaves.size()];
		for(int i=0;i<leaves.size();i++){
			ClusterTree<TALE> leaf = leaves.get( i );
			//System.out.println("adding leaf "+i+" at "+leaf.getOriginalIndex());
			tales[leaf.getOriginalIndex()] = leaf.getClusterElements()[0];
		}
		
		double[][] dmat = computeDistMatrix( tales, costs, at, extraGapOpening, extraGapExtension );
		
		Hclust<TALE> clust = new Hclust<TALE>( null, linkage );
		
		ClusterTree<Integer> superCluster = clust.cluster( dmat, intTrees, -minIndex );
		
		
		ClusterTree<TALEFamily> famTree = superCluster.dropBelow(familyRootOriginalIndexes, families);
		
		return famTree;
		
	}
	
	
	/*private static void leafOrder(ClusterTree<TALE> tree, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		if(tree.getNumberOfElements() == 1 || tree.getSubTrees().length == 1){
			return;
		}else{
			if(tree.getSubTrees().length > 2){
				throw new RuntimeException( "Leaf ordering only implemented for binary ClusterTrees" );
			}else{
				ClusterTree<TALE>[] subs = tree.getSubTrees();
				leafOrder(subs[0],costs,at,extraGapOpening,extraGapExtension);
				leafOrder(subs[1],costs,at,extraGapOpening,extraGapExtension);
				
				TALE t11 = subs[0].getClusterElements()[0];
				TALE t12 = subs[0].getClusterElements()[subs[0].getNumberOfElements()-1];
				TALE t21 = subs[1].getClusterElements()[0];
				TALE t22 = subs[1].getClusterElements()[subs[1].getNumberOfElements()-1];
				
				double d11 = TALEAligner.align( t11, t21, costs, at, extraGapOpening, extraGapExtension ).getCost();
				double d12 = TALEAligner.align( t11, t22, costs, at, extraGapOpening, extraGapExtension ).getCost();
				double d21 = TALEAligner.align( t12, t21, costs, at, extraGapOpening, extraGapExtension ).getCost();
				double d22 = TALEAligner.align( t12, t22, costs, at, extraGapOpening, extraGapExtension ).getCost();
				
				if(d11 < d12 && d11 < d21 && d11 < d22){
					subs[0].reverseOrder();
				}else if(d12 < d21 && d12 < d22){
					subs[0].reverseOrder();
					subs[1].reverseOrder();
				}else if(d22 < d21){
					subs[1].reverseOrder();
				}
				tree.setElements();
			}
		}
		
	}*/
	
	
	private static double[][] computeDistMatrix(TALE[] tales, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		double[][] dmat = new double[tales.length][tales.length];
		
		for(int j=0;j<tales.length;j++){
			for(int k=0;k<tales.length;k++){
				double sc = TALEAligner.align( tales[j], tales[k], costs, at, extraGapOpening, extraGapExtension ).getCost();
				dmat[j][k] = sc;
			}
			//System.out.println(tales[j].getId()+" "+Arrays.toString( dmat[j] ));
		}
		return dmat;
	}
	
	private static double[][] computeDistMatrix2(TALE[] tales, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		double[][] dmat = new double[tales.length][tales.length];
		
		for(int j=0;j<tales.length;j++){
			for(int k=j;k<tales.length;k++){
				double sc = TALEAligner.align( tales[j], tales[k], costs, at, extraGapOpening, extraGapExtension ).getCost();
				dmat[j][k] = sc;
				dmat[k][j] = sc;
			}
			//System.out.println(tales[j].getId()+" "+Arrays.toString( dmat[j] ));
		}
		return dmat;
	}
	
	
	
	private double[][] computeDistMatrix3(TALE[] tales, int firstNew, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		double[][] dmat = new double[tales.length][tales.length];
		
		for(int j=0;j<tales.length;j++){
			for(int k=0;k<tales.length;k++){
				if(j<firstNew && k < firstNew) {
					dmat[j][k] = this.dmat[j][k];
					//dmat[k][j] = this.dmat[k][j];
				}else {
					double sc = TALEAligner.align( tales[j], tales[k], costs, at, extraGapOpening, extraGapExtension ).getCost();
					dmat[j][k] = sc;
					//dmat[k][j] = sc;
				}
			}
			//System.out.println(tales[j].getId()+" "+Arrays.toString( dmat[j] ));
		}
		return dmat;
	}
	
	
	
	
	
	private static Pair<double[][],ClusterTree<TALE>> cluster(TALE[] tales, Linkage linkage, Costs costs, AlignmentType at, double extraGapOpening, double extraGapExtension){
		Hclust<TALE> hclust = new Hclust<TALE>( null, linkage );
		
		double[][] dmat = computeDistMatrix( tales, costs, at, extraGapOpening, extraGapExtension );
		
		//System.out.println();
		
		ClusterTree<TALE> tree = hclust.cluster( dmat, tales );
		
		/*dmat = */
		//tree.leafOrder( dmat );
		/*
		TALE[] els = tree.getClusterElements();
		for(int i=0;i<dmat.length;i++){
			System.out.println(els[i].getId()+" "+Arrays.toString( dmat[i] ));
		}
		System.out.println();
		*/
		
		//leafOrder( tree, costs, at, extraGapOpening, extraGapExtension );
		
		/*ClusterTree<TALE>[] els = tree.getLeaves();
		
		double sum = 0.0;
		for(int i=0;i<els.length-1;i++){
			sum += dmat[ els[i].getOriginalIndex() ][ els[i+1].getOriginalIndex() ];
		}
		
		System.out.println("Sum: "+sum);*/
		
		//tree.leafOrder( dmat );
		
		return new Pair<double[][],ClusterTree<TALE>>(dmat, tree);
	}
	
	
	public void splitClass(String className){
		LinkedList<TALEFamily> famList = new LinkedList<TALEFamilyBuilder.TALEFamily>();

		String[] allFamIDs = ClassAssignmentTool.SchemaFamilyIdGenerator.getFamIDs();
		LinkedList<String> free = new LinkedList<>();
		for(int i=0;i<allFamIDs.length;i++){
			boolean found = false;
			for(int j=0;j<families.length;j++){
				if(allFamIDs[i].equals(families[j].getFamilyId())){
					found = true;
					break;
				}
			}
			if(!found){
				free.add(allFamIDs[i]);
			}
		}
		
		
		for(int i=0;i<families.length;i++){
			if(className.equals(families[i].getFamilyId())){
				double dist = families[i].getTree().getDistance();
				ClusterTree<TALE>[] trees = Hclust.cutTree(dist-1E-6,families[i].getTree());
				trees[0].leafOrder(dmat);
				families[i] = new TALEFamily(families[i].getFamilyId(), trees[0], this);
				renameTALEs(trees[0],families[i].getFamilyId());
				famList.add(families[i]);
				for(int j=1;j<trees.length;j++){
					trees[j].leafOrder(dmat);
					String newName = free.removeFirst();
					renameTALEs(trees[j],newName);
					famList.add(new TALEFamily(newName, trees[j], this));
				}
			}else{
				famList.add(families[i]);
			}
		}
		
		this.families = famList.toArray(new TALEFamily[0]);
		
		LinkedList<TALE> allTALEs = new LinkedList<>();
		for(TALEFamily fam : famList){
			TALE[] tales = fam.getFamilyMembers();
			for(int i=0;i<tales.length;i++){
				allTALEs.add(tales[i]);
			}
		}
		
		dmat = computeDistMatrix(allTALEs.toArray(new TALE[0]), costs, at, extraGapOpening, extraGapExtension);
		
	}
	
	
	private static void renameTALEs(ClusterTree<TALE> clusterTree, String newName) {
		TALE[] tales = clusterTree.getClusterElements();
		ComparableElement<TALE, Integer>[] els = new ComparableElement[tales.length];
		for(int i=0;i<tales.length;i++){
			String id = tales[i].getId();
			id = id.replaceAll(" .*", "");
			id = id.replaceAll("^Tal[A-Z]{2}", "");
			int num = Integer.parseInt(id);
			els[i] = new ComparableElement<TALE, Integer>(tales[i], num);
		}
		Arrays.sort(els);
		for(int i=1;i<=els.length;i++){
			TALE curr = els[i-1].getElement();
			curr.setId( curr.getId().replaceAll("^Tal[A-Z]{2}[0-9]+", newName+i) );
		}
	}

	public void removeTALEsFromFamilies(LinkedList<TALE> toRemove) {
		LinkedList<ClusterTree<TALE>> remainingTALEs = new LinkedList<ClusterTree<TALE>>();
		LinkedList<ClusterTree<TALE>>[] lists = new LinkedList[families.length];
		int minIndex = 0;
		int[] indexMap = new int[getNumberOfTales()];
		for(int i=0;i<families.length;i++){
			ClusterTree<TALE>[] leaves = families[i].getTree().getLeaves();
			for(int j=0;j<leaves.length;j++){
				TALE temp = leaves[j].getClusterElements()[0];
				if(!toRemove.contains(temp)){
					remainingTALEs.add(leaves[j]);
					int idx = leaves[j].getOriginalIndex(); 
					indexMap[idx] = idx;
				}else{
					if(lists[i] == null){
						lists[i] = new LinkedList<ClusterTree<TALE>>();
					}
					int idx = leaves[j].getOriginalIndex();
					indexMap[idx] = -1;
					lists[i].add(leaves[j]);
				}
			}
			if(families[i].getFamilySize() == 0){
				families[i] = null;
			}
			
			int index = families[i].getTree().getMinimumOriginalIndex();
			if(index < minIndex){
				minIndex = index;
			}
		}
		
		for(int i=0,j=0;j<indexMap.length;j++){
			if(indexMap[j] > -1){
				indexMap[j] = i;
				i++;
			}
		}
		
		TALE[] remain = new TALE[remainingTALEs.size()];
		for(int i=0;i<remainingTALEs.size();i++){
			ClusterTree<TALE> rem = remainingTALEs.get(i);
			int idx = rem.getOriginalIndex();
			rem.setOriginalIndex(indexMap[idx]);
			remain[indexMap[idx]] = rem.getClusterElements()[0];
		}
		
		
		double[][] newDmat = computeDistMatrix( remain, costs, at, extraGapOpening, extraGapExtension );//TODO efficiency
		
		LinkedList<TALEFamily> famList = new LinkedList<TALEFamilyBuilder.TALEFamily>();
		
		for(int i=0;i<lists.length;i++){
			if(lists[i] != null){
				TALEFamily temp = families[i].removeTALE(newDmat, lists[i], -minIndex, this);
				if(temp != null){
					famList.add(temp);
				}
			}else{
				famList.add(families[i]);
			}
		}
		
		this.families = famList.toArray(new TALEFamily[0]);
		this.dmat = newDmat;
		
	}
	
	private int getNumberOfTales() {
		int num = 0;
		for(int i=0;i<families.length;i++){
			num += families[i].getFamilySize();
		}
		return num;
	}

	public void addTALEsToFamilies(Pair<Integer,LinkedList<TALE>>[] assignment, TALE[] unassigned, FamilyIdGenerator idGenerator){
		
		if(idGenerator == null){
			idGenerator = new DefaultFamilyIdGenerator();
		}
		
		int numNew = unassigned.length;
		for(int i=0;i<assignment.length;i++){
			numNew += assignment[i].getSecondElement().size();
		}
		
		TALE[] allTALEs = new TALE[dmat.length + numNew];
		
		int minIndex = 0;
		int maxIndex = 0;
		for(int i=0;i<families.length;i++){
			ClusterTree<TALE> famTree = families[i].getTree();
			
			ClusterTree<TALE>[] leaves = famTree.getLeaves();
			for(int j=0;j<leaves.length;j++){
				int idx = leaves[j].getOriginalIndex(); 
				allTALEs[idx] = leaves[j].getClusterElements()[0];
				if(idx > maxIndex){
					maxIndex = idx;
				}
			}
			
			int index = famTree.getMinimumOriginalIndex();
			if(index < minIndex){
				minIndex = index;
			}
			
		}
		
		if(maxIndex+1 != dmat.length){
			throw new RuntimeException("Indexes do not match "+(maxIndex+1)+" <-> "+dmat.length);
		}
		
		int k=maxIndex+1;
		for(int i=0;i<assignment.length;i++){
			LinkedList<TALE> nt = assignment[i].getSecondElement();
			for(int j=0;j<nt.size();j++,k++){
				allTALEs[k] = nt.get( j );
			}
		}
		for(int i=0;i<unassigned.length;i++,k++){
			allTALEs[k] = unassigned[i];
		}
		
		//double[][] newDmat2 = computeDistMatrix( allTALEs, costs, at, extraGapOpening, extraGapExtension );//TODO efficiency
		double[][] newDmat = computeDistMatrix3( allTALEs, maxIndex+1, costs, at, extraGapOpening, extraGapExtension );//TODO efficiency
		
		/*
		 * for(int i=0;i<newDmat.length;i++) { for(int j=0;j<newDmat[i].length;j++) {
		 * if(Math.abs(newDmat[i][j] - newDmat2[i][j])>1E-10) {
		 * System.out.print(i+","+j+": "+newDmat[i][j]+" <-> "+newDmat2[i][j]);
		 * if(i<this.dmat.length && j<this.dmat[i].length) {
		 * System.out.print(" <- "+this.dmat[i][j]); } System.out.println(
		 * " "+allTALEs[i].getId()+" "+ allTALEs[j].getId()); } } }
		 */
		
		TALEFamily[] newFams = new TALEFamily[families.length];
		System.arraycopy( families, 0, newFams, 0, families.length );
		
		k=maxIndex+1;
		for(int i=0;i<assignment.length;i++){
			Integer famIdx = assignment[i].getFirstElement();
			LinkedList<TALE> newTALEs = assignment[i].getSecondElement();
			for(int j=0;j<newTALEs.size();j++,k++){
				newFams[ famIdx ] = newFams[ famIdx ].addTALE( newDmat, newTALEs.get(j), k, -minIndex, this );
				int newMinIdx = newFams[ famIdx ].tree.getMinimumOriginalIndex();
				if(newMinIdx < minIndex){
					minIndex = newMinIdx;
				}
			}
		}

		ClusterTree<TALE>[] createdTrees = null;
		
		if(unassigned.length > 0){

			ClusterTree<TALE>[] newLeaves = new ClusterTree[unassigned.length];
			for(int i=0;i<unassigned.length;i++,k++){
				newLeaves[i] = new ClusterTree<TALE>(unassigned[i],k);
			}

			Hclust<TALE> hclust = new Hclust<TALE>( null, linkage );

			ClusterTree<TALE> newTree = hclust.cluster( -minIndex, newDmat, newLeaves );

			createdTrees = Hclust.cutTree( cut, newTree );

		}else{
			createdTrees = new ClusterTree[0];
		}
		
		TALEFamily[] allFams = new TALEFamily[newFams.length + createdTrees.length];
		System.arraycopy( newFams, 0, allFams, 0, newFams.length );
		for(int i=0;i<createdTrees.length;i++){
			createdTrees[i].leafOrder(newDmat);
			allFams[ i + newFams.length ] = new TALEFamily( null, createdTrees[i], this );
		}
		
		idGenerator.setFamilyIDs( allFams, this );
		
		this.dmat = newDmat;
		this.families = allFams;
		

	}
	
	
	public Pair<TALEFamily,Double> getClosestFamily(TALE tale, FamilyDistance dist ){
		Pair<Integer,Double> pair = getClosestFamilyIndex( tale, dist );
		return new Pair<TALEFamily,Double>(getFamily( pair.getFirstElement() ),pair.getSecondElement());
	}
	
	public Pair<Integer,Double> getClosestFamilyIndex(TALE tale, FamilyDistance dist){
		int closest = -1;
		double min = Double.POSITIVE_INFINITY;
		for(int i=0;i<this.families.length;i++){
			double d = this.families[i].getDistance( tale, dist, this );
			if(d < min){
				min = d;
				closest = i;
			}
		}
		/*int n = this.families[closest].getFamilySize();//TODO start
		if(this.linkage == Linkage.AVERAGE){
			min = ((n*(n-1))/2* this.families[closest].getTree().getDistance() + n*min)/((n*(n+1))/2);
		}else if(this.linkage == Linkage.COMPLETE){
			min = Math.max(min, this.families[closest].getTree().getDistance());
		}*///TODO end
		return new Pair<Integer,Double>(closest,min);
	}
	
	public Pair<TALEFamily,Double> getMostSignificantFamily(TALE tale, AlignmentPValues pv, FamilyDistance dist){
		Pair<Integer,Double> pair = getMostSignificantFamilyIndex( tale, pv, dist );
		return new Pair<TALEFamily,Double>( getFamily( pair.getFirstElement()), pair.getSecondElement() );
	}
	
	public Pair<Integer,Double> getMostSignificantFamilyIndex(TALE tale, AlignmentPValues pv, FamilyDistance dist){
		int closest = -1;
		double min = Double.POSITIVE_INFINITY;
		for(int i=0;i<this.families.length;i++){
			double d = this.families[i].getSignificance( tale, pv, dist, this );
			if(d < min){
				min = d;
				closest = i;
			}
		}
		return new Pair<Integer,Double>(closest,min);
	}
	
	public TALEFamily getFamily(int index){
		return this.families[index];
	}
	
	public TALEFamily[] getFamilies(){
		return this.families;
	}

	
	public String[] getReservedNames() {
		return reservedNames;
	}

	
	public void setReservedNames( String[] reservedNames ) {
		this.reservedNames = reservedNames;
	}

	public Costs getCosts() {
		return costs;
	}

	public double getPVal() {
		return pval;
	}

	
	
}
