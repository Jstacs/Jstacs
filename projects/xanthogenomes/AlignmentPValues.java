package projects.xanthogenomes;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Set;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.graphics.PDFAdaptor;
import projects.xanthogenomes.BuildFamilies.FamilyResult;
import projects.xanthogenomes.Tools.Translator;
import projects.xanthogenomes.alignmentCosts.RVDCosts;


public class AlignmentPValues {

	private double[][] cost;
	private double[][] prob;
	
	private static int n=0;
	
	private static double doublePrec = 1E-3;
	
	public AlignmentPValues(TALE[] allTales, RVDCosts costs){
		
		AlphabetContainer rvdAlph = allTales[0].getRvdSequence().getAlphabetContainer();
		
		double[] rvdProbs = new double[(int)rvdAlph.getAlphabetLengthAt( 0 )];
		Arrays.fill( rvdProbs, 1.0/rvdProbs.length );
		
		for(int i=0;i<allTales.length;i++){
			Sequence rvds = allTales[i].getRvdSequence();
			for(int j=0;j<rvds.getLength();j++){
				rvdProbs[ rvds.discreteVal( j ) ] ++;
			}
		}
		
		Normalisation.sumNormalisation( rvdProbs );
		
		
		cost = new double[rvdProbs.length][];
		prob = new double[rvdProbs.length][];
		
		for(int i=0;i<rvdProbs.length;i++){
			HashMap<Double,Double> map = new HashMap<Double,Double>();
			for(int j=0;j<rvdProbs.length;j++){
				Double c = costs.getCostFor( rvdAlph, i, j );
				if(map.containsKey( c )){
					map.put( c, map.get( c ) + rvdProbs[j] );
				}else{
					map.put( c, rvdProbs[j] );
				}
			}
			
			Set<Double> keys = map.keySet();
			LinkedList<Double> list = new LinkedList<Double>( keys );
			Collections.sort( list );
			
			cost[i] = new double[list.size()];
			prob[i] = new double[list.size()];
			for(int j=0;j<list.size();j++){
				cost[i][j] = list.get( j );
				prob[i][j] = Math.log( map.get( list.get( j ) ) );
			}
			//System.out.println(i+" "+Arrays.toString( cost[i] ));
			
		}
		
	}	
	
	
	/*public double getLog10PValue( Sequence rvds, double costThresh ){
		int[] costIndex = new int[rvds.getLength()];
		
		double startp = 0;
		double startCost = 0;
		for(int i=0;i<rvds.getLength();i++){
			int rvd = rvds.discreteVal( i );
			startp += prob[rvd][0];
			startCost += cost[rvd][0];
		}
		
		double p2 = Double.NEGATIVE_INFINITY;
		if(startCost <= costThresh){
			p2 = getLogPValue( rvds, costThresh, costIndex, 0, startCost, startp );
		}
		
		
		return Normalisation.getLogSum( p2, startp )/Math.log( 10 );
		
	}
	
	private double getLogPValue( Sequence rvds, double costThresh, int[] costIndex, int off, double currScore, double currProb ){
		
		int cl = cost[0].length;
		int sl = rvds.getLength();
		
		double[] temp = new double[ (sl-off)*cost[0].length*2 + 1 ];
		temp[0] = Double.NEGATIVE_INFINITY;
		int k=1;
		for(int i=sl-1;i>=off;i--){
			
			int rvd = rvds.discreteVal( i );

			double tempScore = currScore - cost[rvd][0];
			double tempProb = currProb - prob[rvd][0];
			
			for(int j=1;j<cl;j++){
				//costIndex[i] = j;
				double myScore = tempScore + cost[rvd][j];
				//System.out.println(i+" "+j+" "+rvd+" "+myScore);
				if(myScore <= costThresh + doublePrec){
					
					n++;
					double myProb = tempProb + prob[rvd][j];
					
					temp[k] = myProb;
					k++;
					//costIndex[i] = j;
					
					//System.out.println(myProb);
					//System.out.println(Arrays.toString( costIndex ));
					double p2 = getLogPValue( rvds, costThresh, costIndex, i+1, myScore, myProb );
					temp[k] = p2;
					k++;					
					
				}
				//System.out.println();
			}
			//costIndex[i] = 0;
		}
		
		double logp = Normalisation.getLogSum( 0, k, temp );
		
		return logp;
	}*/
	
	
	
	public double getLog10PValue( Sequence rvds, double costThresh, double baseScore ){
		
		double[] prevScores = new double[]{baseScore};
		double[] prevProbs = new double[]{0};
				
		for(int i=0;i<rvds.getLength();i++){
			int rvd = rvds.discreteVal( i );
			double[] sc = cost[rvd];
			double[] pr = prob[rvd];
			
			Pair<double[],double[]> pair = joinAggregateAndFilter(prevScores,prevProbs,sc,pr,costThresh);
			
			prevScores = pair.getFirstElement();
			prevProbs = pair.getSecondElement();
			if(prevScores.length == 0){
				return Double.NEGATIVE_INFINITY;
			}
		//	System.out.println(i+" "+prevScores.length);
		//	System.out.println( Arrays.toString(prevScores) );
		//	System.out.println( Arrays.toString(prevProbs) );
		}
		
		int i=0;
		while( i<prevScores.length && prevScores[i] < costThresh ){
			i++;
		}
		
		return Normalisation.getLogSum( 0,i,prevProbs )/Math.log( 10 );
		
		
	}
	
	
	
	private Pair<double[], double[]> joinAggregateAndFilter( double[] prevScores, double[] prevProbs, double[] sc, double[] pr,
			double costThresh ) {
		
		double[] tempSc = new double[prevScores.length*sc.length];
		double[] tempPr = new double[prevProbs.length*pr.length];
		for(int i=0,k=0;i<prevScores.length;i++){
			for(int j=0;j<sc.length;j++,k++){
				tempSc[k] = prevScores[i] + sc[j];
				tempPr[k] = prevProbs[i] + pr[j];
			}
		}
		
		ToolBox.sortAlongWith( tempSc, tempPr );
		
		DoubleList resSc = new DoubleList(tempSc.length);
		DoubleList resPr = new DoubleList(tempSc.length);
		int i=0;
		//System.out.println(tempSc.length);
		while(i<tempSc.length){
			
			int j=i;
			double prev = tempSc[j];
			/*if(prev-costThresh > doublePrec){
				break;
			}*/
			while(i<tempSc.length && Math.abs( tempSc[i] - prev ) < doublePrec){
				i++;
			}
			
			resSc.add( prev );
			//System.out.println(j+" "+i);
			resPr.add( Normalisation.getLogSum( j,i,tempPr ) );
		}
		//System.out.println();
		
		return new Pair<double[],double[]>(resSc.toArray(),resPr.toArray());
	}


	public static void main(String[] args) throws Exception {
		
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder( FileManager.readFile(args[0]) );
		
		builder.setToOld();
		
		TALE[] tales = builder.getAllTALEs();
		
		Pair<TALEFamilyBuilder,FamilyResult[]> fams = BuildFamilies.build(tales, Integer.MAX_VALUE, 1.1);
		
		FamilyResult res = fams.getSecondElement()[0];
		
		TALEFamilyTreePlotter plotter = new TALEFamilyTreePlotter(10);
		
		PDFAdaptor ad = new PDFAdaptor();
		
		int[] dim = plotter.getDimension(ad.getGraphics(10, 10), res.getFamily());
		
		ad = new PDFAdaptor();
		
		plotter.plot(ad.getGraphics(dim[0], dim[1]), res.getFamily());
		
		ad.generateOutput("/Users/dev/Downloads/completeTree.pdf");
		
		
		/*TALE longest = tales[0];
		for(int i=1;i<tales.length;i++){
			if(tales[i].getNumberOfRepeats() > longest.getNumberOfRepeats()){
				longest = tales[i];
			}
		}
		
		RVDCosts costs = new RVDCosts( 1.0, 0.2, 0.8, 0.0 );//TODO
		
		Sequence rvd = longest.getRvdSequence();
		rvd = rvd.getSubSequence(6);
		AlphabetContainer rvdAlph = TALE.RVDAlphabet;
		AlignmentPValues pv = new AlignmentPValues( tales, costs );
		System.out.println(rvd);
		pv.getLog10PValue(rvd, Double.POSITIVE_INFINITY, 0.0);*/
		
	}


	public double getLog10PValue( TALE tale, TALE tale2, double sc, double extraGapOpening, double extraGapExtension ) {
		
		String str = null;
		//tale should always be the shorter one
		if(tale.getNumberOfRepeats() > tale2.getNumberOfRepeats()){
			TALE temp = tale;
			tale = tale2;
			tale2 = temp;
			
		}
		
		Sequence rvd = tale.getRvdSequence();
		Sequence rvd2 = tale2.getRvdSequence();
		
		double sum = 0.0;
		double sumq = 0.0;
		int l = rvd2.getLength()-rvd.getLength()+1;
		for(int i=0;i<l;i++){
			double base = ((i>0 && l>1) ? extraGapOpening : 0) + (i<l-1 && l>1 ? extraGapOpening : 0) + (l-1)*extraGapExtension;
			//System.out.println(rvd+"\n"+rvd2+" "+base+" "+sc);
			double p = this.getLog10PValue( rvd2.getSubSequence( i, rvd.getLength() ), sc, base );//TODO gaps
			p = log1m( p * Math.log( 10 ) );
			sum += p;
			
			double q = this.getLog10PValue( rvd, sc, base );//TODO gaps
			q = log1m(  q * Math.log( 10 ) );
			sumq += q;
		}
		
		
		
		//double p1 = log1m( sum );
		//double p2 = log1m( sumq);
		
		double ret = log1m( sum + sumq )/Math.log( 10 );
		//System.out.println(p1+" "+p2+" "+ret);
		
		return ret;
	}
	
	
	private static double log1m(double d){
		
		double d1 = Math.log1p( -Math.exp(d) );//TODO improve precision
		
		return d1;
	}
	
	
	
}
