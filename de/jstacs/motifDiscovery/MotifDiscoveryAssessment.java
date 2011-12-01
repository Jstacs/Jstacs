/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.motifDiscovery;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.sequences.annotation.LocatedSequenceAnnotation;
import de.jstacs.data.sequences.annotation.LocatedSequenceAnnotationWithLength;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.results.ListResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.utils.DoubleList;


/**
 * This class enables the user to assess the prediction of  motif occurrences
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class MotifDiscoveryAssessment {

	/**
	 * This method computes the nucleotide and site measures.
	 * 
	 * @param truth the {@link Sample} annotated with the true annotation
	 * @param prediction  annotated with the predicted annotation
	 * @param maxDiff the maximal difference between predicted and true start position;
	 * 		this value is used to determine the site measures
	 * 
	 * @return a {@link ListResult} containing all {@link NumericalResultSet}s
	 * 
	 * @throws Exception if something went wrong
	 */
	public static ListResult assess(Sample truth, Sample prediction, int maxDiff) throws Exception{
		if(truth.getNumberOfElements() != prediction.getNumberOfElements()){
			throw new Exception("Truth must contain same number of sequences as prediction.");
		}
		
		HashMap<String,IdentifiersOfType> mapT = new HashMap<String, IdentifiersOfType>();
		HashMap<String,IdentifiersOfType> mapP = new HashMap<String, IdentifiersOfType>();
		
		Sequence seqT, seqP;
		for(int i=0;i<truth.getNumberOfElements();i++){
			seqT = truth.getElementAt( i );
			seqP = prediction.getElementAt( i );
			if( !seqT.equals( seqP ) ) {
				throw new Exception( "samples incomparable: check sequence " + i );
			}
			
			SequenceAnnotation[] ann = seqT.getAnnotation();
			if(ann != null){
				fillLocatedSequenceAnnotations( mapT, ann );
			}
			ann = seqP.getAnnotation();
			if(ann != null){
				fillLocatedSequenceAnnotations( mapP, ann );
			}
		}
		
		NumericalResultSet[] ress = new NumericalResultSet[mapT.size()];
		
		Collection<IdentifiersOfType> collT = mapT.values();
		Iterator<IdentifiersOfType> itT = collT.iterator();
		
		int idx = 0;
		while(itT.hasNext()){
			IdentifiersOfType idT = itT.next();
			idT.reset();
			idT.nextPermutation();
			int[][] annotT = new int[truth.getNumberOfElements()][];
			for(int i=0;i<truth.getNumberOfElements();i++){
				annotT[i] = getAnnot( truth.getElementAt( i ), idT );
			}
			IdentifiersOfType idP = mapP.get( idT.getType() );
			if(idP == null){
				ress[idx] = new NumericalResultSet();
				idx++;
				continue;
			}
			idP.reset();
			int[][] annotP = new int[prediction.getNumberOfElements()][];
			
			NumericalResultSet bestNMeasures = null;
			NumericalResultSet bestSMeasures = null;
			
			while(idP.nextPermutation()){
				for(int i=0;i<prediction.getNumberOfElements();i++){
					annotP[i] = getAnnot( prediction.getElementAt( i ), idP );
				}
				
				NumericalResultSet nMeasures = computeNucleotideLevelMeasures( annotT, annotP, idT );
				if(bestNMeasures == null || bestNMeasures.getResultAt( bestNMeasures.getNumberOfResults() - 1 ).getResult().compareTo( nMeasures.getResultAt( nMeasures.getNumberOfResults() - 1 ).getResult() ) < 0){
					bestNMeasures = nMeasures;
					bestSMeasures = computeSiteLevelMeasures( truth, prediction, maxDiff, idT, idP );
				}	
			}
			ress[idx] = new NumericalResultSet( ArrayHandler.cast( NumericalResult.class, bestNMeasures.getResults() ),
					ArrayHandler.cast( NumericalResult.class, bestSMeasures.getResults()));
			idx++;
		}
		return new ListResult("Motif discovery assessment","Assessment of the motif discovery result on nucleotide and site level",null,ress);
	}
	
	private static NumericalResultSet computeSiteLevelMeasures(Sample truth, Sample prediction, int maxDiff, IdentifiersOfType idT, IdentifiersOfType idP){
		LinkedList<LocatedSequenceAnnotation> listT = new LinkedList<LocatedSequenceAnnotation>();
		LinkedList<LocatedSequenceAnnotation> listP = new LinkedList<LocatedSequenceAnnotation>();
		
		LinkedList<NumericalResult> list = new LinkedList<NumericalResult>();
		
		Iterator<String> it = idT.keySet().iterator();
		while(it.hasNext()){
			double tp=0,fn=0,fp=0;
			String t = it.next();
			for(int i=0;i<truth.getNumberOfElements();i++){
				listT.clear();
				listP.clear();
				Sequence seq = truth.getElementAt( i );
				getLSA( listT, seq.getAnnotation(), idT.getType(), t );
				
				seq = prediction.getElementAt( i );
				getLSA( listP, seq.getAnnotation(), idP.getType(), idP.getIdentifierForIdx( idT.get( t ) ) );
				
				Iterator<LocatedSequenceAnnotation> itT = listT.iterator();
				Iterator<LocatedSequenceAnnotation> itP = null;
				while(itT.hasNext()){
					boolean found = false;
					LocatedSequenceAnnotation anT = itT.next();
					itP = listP.iterator();
					while(itP.hasNext()){
						LocatedSequenceAnnotation anP = itP.next();
						if(Math.abs(anT.getPosition() - anP.getPosition()) <= maxDiff){
							found = true;
						}
					}
					if(found){
						tp++;
					}else{
						fn++;
					}
				}
				
				itP = listP.iterator();
				while(itP.hasNext()){
					boolean found = false;
					LocatedSequenceAnnotation anP = itP.next();
					itT = listT.iterator();
					while(itT.hasNext()){
						LocatedSequenceAnnotation anT = itT.next();
						if(Math.abs(anT.getPosition() - anP.getPosition()) <= maxDiff){
							found = true;
						}
					}
					if(!found){
						fp++;
					}
				}
				
			}
			list.add( new NumericalResult("sSn: "+idT.getType()+", "+t,"Sensitivity for type "+idT.getType()+" and identifier "+t+" (site level)",(double)tp/(double)(tp+fn)) );
			list.add( new NumericalResult("sPPV: "+idT.getType()+", "+t,"Positive predictive value for type "+idT.getType()+" and identifier "+t+" (site level)",(double)tp/(double)(tp+fp)) );
			list.add( new NumericalResult("sPC: "+idT.getType()+", "+t,"Performance coefficient for type "+idT.getType()+" and identifier "+t+" (site level)",(double)tp/(double)(tp+fn+fp)) );
			list.add( new NumericalResult("sASP: "+idT.getType()+", "+t,"Average site performance for type "+idT.getType()+" and identifier "+t+" (site level)",
					( ((double)tp/(double)(tp+fn)) + ((double)tp/(double)(tp+fp)) )/2.0) );
		}
		return new NumericalResultSet(list);
	}
	
	private static void getLSA(LinkedList<LocatedSequenceAnnotation> list, SequenceAnnotation[] anns, String type, String id){
		if(anns == null){
			return;
		}
		for(int i=0;i<anns.length;i++){
			SequenceAnnotation ann = anns[i];
			if(ann instanceof LocatedSequenceAnnotation && ann.getType().equals( type ) && ann.getIdentifier().equals( id )){
				list.add( (LocatedSequenceAnnotation) ann);
			}
			SequenceAnnotation[] sub = ann.getSubAnnotations();
			if(sub != null){
				getLSA( list, sub, type, id );
			}
		}
	}
	
	private static NumericalResultSet computeNucleotideLevelMeasures(int[][] truth, int[][] prediction, IdentifiersOfType ids){
		Iterator<Map.Entry<String, Integer>> it = ids.entrySet().iterator();
		LinkedList<NumericalResult> list = new LinkedList<NumericalResult>();
		while(it.hasNext()){
			Map.Entry<String, Integer> en = it.next();
			double tp=0,tn=0,fp=0,fn=0;
			for(int i=0;i<truth.length;i++){
				for(int j=0;j<truth[i].length;j++){
					if(truth[i][j] == en.getValue()){
						if(prediction[i][j] == truth[i][j]){
							tp ++;
						}else{
							fn ++;
						}
					}else{
						if(prediction[i][j] != en.getValue()){
							tn ++;
						}else{
							fp ++;
						}
					}
				}
			}
			list.add( new NumericalResult("nSn: "+ids.getType()+", "+en.getKey(),"Sensitivity for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",(double)tp/(double)(tp+fn)) );
			list.add( new NumericalResult("nSp: "+ids.getType()+", "+en.getKey(),"Specificity for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",(double)tn/(double)(tn+fp)) );
			list.add( new NumericalResult("nPPV: "+ids.getType()+", "+en.getKey(),"Positive predictive value for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",(double)tp/(double)(tp+fp)) );
			list.add( new NumericalResult("nPC: "+ids.getType()+", "+en.getKey(),"Performance coefficient for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",(double)tp/(double)(tp+fn+fp)) );
			list.add( new NumericalResult("nASP: "+ids.getType()+", "+en.getKey(),"Average site performance for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",
					( ((double)tp/(double)(tp+fn)) + ((double)tp/(double)(tp+fp)) )/2.0) );
			list.add( new NumericalResult("nCC: "+ids.getType()+", "+en.getKey(),"Correlation coefficient for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",
					(tp*tn - fn*fp) / Math.sqrt( (tp + fn) * (tn + fp) * (tp + fp) * (tn + fn) ) ) );
			list.add( new NumericalResult("nCR: "+ids.getType()+", "+en.getKey(),"Classification rate for type "+ids.getType()+" and identifier "+en.getKey()+" (nucleotide level)",
					(double)(tp+tn)/(double)(tp+tn+fp+fn)) );			
		}
		int corr=0,wr=0;
		for(int i=0;i<truth.length;i++){
			for(int j=0;j<truth[i].length;j++){
				if(truth[i][j] == prediction[i][j]){
					corr ++;
				}else{
					wr ++;
				}
			}
		}
		list.add( new NumericalResult("ngCR: "+ids.getType(),"Global classification rate for type "+ids.getType()+" (nucleotide level)",(double)corr/(double)(corr+wr)) );
		return new NumericalResultSet(list);
	}
	
	private static int[] getAnnot(Sequence seq, IdentifiersOfType ids){
		int[] annot = new int[seq.getLength()];
		Arrays.fill( annot, -1 );
		SequenceAnnotation[] ann = seq.getAnnotation();
		if(ann != null){
			getAnnot( annot, ann, ids );
		}
		return annot;
	}
	
	private static void getAnnot(int[] annot, SequenceAnnotation[] anns, IdentifiersOfType ids){
		for(int i=0;i<anns.length;i++){
			SequenceAnnotation ann = anns[i];
			if(ann.getType().equals( ids.getType() )){
				if(ann instanceof LocatedSequenceAnnotation){
					int pos = ((LocatedSequenceAnnotation) ann).getPosition();

					int idx = ids.getIdxForIdentifier( (LocatedSequenceAnnotation) ann );
					annot[pos] = idx;
					if(ann instanceof LocatedSequenceAnnotationWithLength){
						int length = ((LocatedSequenceAnnotationWithLength) ann).getLength();
						for(int j=1;j<length;j++){
							annot[pos+j] = idx;
						}
					}
				}
			}
			SequenceAnnotation[] sub = ann.getSubAnnotations();
			if(sub != null){
				getAnnot( annot, sub, ids );
			}
		}
	}
	
	private static void fillLocatedSequenceAnnotations(HashMap<String,IdentifiersOfType> map, SequenceAnnotation[] anns ) throws Exception{
		for(int i=0;i<anns.length;i++){
			SequenceAnnotation ann = anns[i];
			if(ann instanceof LocatedSequenceAnnotation){
				if(map.get( ann.getType() ) == null){
					map.put( ann.getType(), new IdentifiersOfType(ann.getType()) );
				}
				map.get( ann.getType() ).addIdentifier( (LocatedSequenceAnnotation) ann );
			}
			SequenceAnnotation[] sub = ann.getSubAnnotations();
			if( sub != null){
				fillLocatedSequenceAnnotations( map, sub );
			}
		}
	}
	
	private static class IdentifiersOfType extends HashMap<String, Integer>{
		
		private String type;
		private int nextIdx;
		private int[] perm;
		private int exchangeNext;
		
		public IdentifiersOfType(String type){
			super();
			this.type = type;
			this.nextIdx = 0;
		}
		
		public String getType(){
			return type;
		}
		
		public void addIdentifier(LocatedSequenceAnnotation ann) throws Exception{
			if(ann.getType() != type){
				throw new Exception("Wrong type");
			}
			String id = ann.getIdentifier();
			if(this.get( id ) == null){
				this.put( id, nextIdx );
				nextIdx++;
			}
		}
		
		public String getIdentifierForIdx(int idx){
			Iterator<Map.Entry<String, Integer>> it = this.entrySet().iterator();
			while(it.hasNext()){
				Map.Entry<String, Integer> en = it.next();
				if(en.getValue() == idx){
					return en.getKey();
				}
			}
			return null;
		}
		
		public int getIdxForIdentifier(LocatedSequenceAnnotation ann){
			String id = ann.getIdentifier();
			return perm[this.get( id )];
		}
		
		public void reset(){
			perm = null;
		}
		
		public boolean nextPermutation(){
			if(perm == null){
				perm = new int[this.size()];
				for(int i=0;i<perm.length;i++){
					perm[i] = i;
				}
				exchangeNext = perm.length - 1;
				return true;
			}else if(perm.length == 1){
				return false;
			}else{
				int temp = perm[exchangeNext];
				perm[exchangeNext] = perm[exchangeNext - 1];
				perm[exchangeNext - 1] = temp;
				
				exchangeNext --;
				if(exchangeNext == 0){
					exchangeNext = perm.length - 1;
				}
				
				for(int i=0;i<perm.length-1;i++){
					if(perm[i] != perm[i+1]-1){
						return true;
					}
				}
				return false;
			}
		}
		
	}
	
	/**
	 * This method provides some score arrays that can be used in {@link de.jstacs.classifier.performanceMeasures.AbstractPerformanceMeasure} to determine some
	 * curves or area under curves based on the values of the predictions. The scores are generated by <code>offset+factor*values[i][j]</code>.
	 * 
	 * @param data the data
	 * @param values the array of smoothed values
	 * @param offset the offset that is added to the current values
	 * @param factor the factor that is multiplied to the current values 
	 * @param identifier the identifier of the annotation of the positive class
	 * 
	 * @return two arrays containing values; the first contains the values for the positive class, the second contains
	 * 			the values for the negative class
	 * 
	 * @see SequenceAnnotation#getIdentifier()
	 */
	public static double[][] getSortedValuesForMotifAndFlanking( Sample data, double[][] values, double offset, double factor, String identifier ) {
		boolean[] added = null;
		DoubleList pos = new DoubleList(1000);
		DoubleList neg = new DoubleList(100000);
		Sequence seq;
		SequenceAnnotation[] seqAn;
		LocatedSequenceAnnotationWithLength lsawl;
		for( int s, l, k, j, i = 0; i < data.getNumberOfElements(); i++ ) {
			seq = data.getElementAt(i);
			seqAn = seq.getAnnotation();
			if( added == null || added.length != seq.getLength() ) {
				added = new boolean[seq.getLength()];
			}
			Arrays.fill( added, false );
			if( seqAn != null ) {
				for( j = 0; j < seqAn.length; j++ ) {
					if( seqAn[j] instanceof LocatedSequenceAnnotationWithLength
							&& seqAn[j].getIdentifier().equalsIgnoreCase( identifier ) ) {
						lsawl = (LocatedSequenceAnnotationWithLength) seqAn[j];
						s = lsawl.getPosition();
						l = lsawl.getLength();
						for( k = 0; k < l; k++, s++ ) {
							if( !added[s] ){
								pos.add( offset + factor*values[i][s] );
								added[s] = true;
							}
						}
					}
				}
			}
			for( k = 0; k < added.length; k++ ) {
				if( !added[k] ) {
					neg.add( offset + factor*values[i][k] );
				}
			}
		}
		double[][] res = { pos.toArray(), neg.toArray() };
		Arrays.sort( res[0] );
		Arrays.sort( res[1] );
		return res;
	}
	
	
	/**
	 * Returns the scores read from the prediction <code>pred</code> for the motif with identifier <code>identifier</code> and flanking sequences as annotated in
	 * the {@link Sample} data. The <code>identifier</code> may be <code>null</code> to obtain the scores for all motifs, irrespective of present identifiers.
	 * The first dimension of the returned array contains the scores for the motif annotations, while the second dimension contains the scores of the flanking sequences. 
	 * Both dimensions are sorted and can be directly used in the methods of {@link de.jstacs.classifier.performanceMeasures.AbstractPerformanceMeasure}.
	 * The scores for the predictions must be added to the {@link LocatedSequenceAnnotationWithLength} representing the motifs as additional annotation using {@link LocatedSequenceAnnotationWithLength#LocatedSequenceAnnotationWithLength(String, String, LocatedSequenceAnnotation[], Result...)}
	 * with the name of the annotation, i.e. the name of the corresponding {@link Result} equal to &quot;score&quot;.
	 * @param data the {@link Sample} annotated with the truth
	 * @param pred the {@link Sample} annotated with the prediction and associated scores
	 * @param identifier the identifier of the motif
	 * @return the scores of motifs and flanking sequences
	 */
	public static double[][] getSortedScoresForMotifAndFlanking( Sample data, Sample pred, String identifier ) {
		boolean[] added = new boolean[data.getElementLength()];
		double[] scores = new double[data.getElementLength()];
		DoubleList pos = new DoubleList(1000);
		DoubleList neg = new DoubleList(100000);
		SequenceAnnotation[] seqAn;
		SequenceAnnotation[] predAn;
		Result[] addRes;
		LocatedSequenceAnnotationWithLength lsawl;
		for( int s, l, k, j, i = 0; i < data.getNumberOfElements(); i++ ) {
			Arrays.fill( added, false );
			seqAn = data.getElementAt(i).getAnnotation();
			predAn = pred.getElementAt( i ).getAnnotation();
			Arrays.fill( scores, -Double.MAX_VALUE );
			for(int p=0;predAn != null && p<predAn.length;p++){
				if(predAn[p] instanceof LocatedSequenceAnnotationWithLength){
					lsawl = (LocatedSequenceAnnotationWithLength) predAn[p];
					int start = lsawl.getPosition();
					int end = lsawl.getEnd();
					double score = -Double.MAX_VALUE;
					addRes = lsawl.getAnnotations();
					for(int q=0;q<addRes.length;q++){
						if("score".equals( addRes[q].getName()) && addRes[q] instanceof NumericalResult ){
							score = (Double) ((NumericalResult)addRes[q]).getResult();
						}
					}
					for(int q=start;q<end;q++){
						scores[q] = Math.max( score, scores[q]);
					}
				}
			}
			
			if( seqAn != null ) {
				for( j = 0; j < seqAn.length; j++ ) {
					if( seqAn[j] instanceof LocatedSequenceAnnotationWithLength
							&& seqAn[j].getIdentifier().equalsIgnoreCase( identifier ) ) {
						lsawl = (LocatedSequenceAnnotationWithLength) seqAn[j];
						s = lsawl.getPosition();
						l = lsawl.getLength();
						for( k = 0; k < l; k++, s++ ) {
							if( !added[s] ){
								pos.add( scores[s] );
								added[s] = true;
							}
						}
					}
				}
			}
			for( k = 0; k < added.length; k++ ) {
				if( !added[k] ) {
					neg.add( scores[k] );
				}
			}
		}
		double[][] res = { pos.toArray(), neg.toArray() };
		Arrays.sort( res[0] );
		Arrays.sort( res[1] );
		return res;
	}
	
	
}
