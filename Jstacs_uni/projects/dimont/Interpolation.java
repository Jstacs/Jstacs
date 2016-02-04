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

package projects.dimont;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.WeightedDataSetFactory;
import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.ToolBox.TiedRanks;

public enum Interpolation {
	LOG_LINEAR,
	LINEAR_jens,
	LINEAR,
	RANK_OUTLIER,
	RANK_OUTLIER_REPEATED,
	RANK_LOG,
	RANK_SIGMOID,
	SIGMOID,
	PERCENTILE_LOGISTIC,
	RANK_SELEX;
	
	public static double[] getWeight( DataSet data, double[] intensa, double aPrioriWeight, Interpolation interpol ) throws Exception {
		double[] weights = new double[intensa.length];
		
		double min, max, thresh;
		double[] help = intensa.clone();
		Arrays.sort( help );
		min = help[0];
		max = help[help.length-1];
		thresh = help[ (int) Math.ceil( (1.0-aPrioriWeight)*help.length )];
		
		double a1, b1, c = 1E-10, h, m=0, q25, q75, temp;
		int idx = -1;
		int[] rank;
		switch(interpol) {
			case PERCENTILE_LOGISTIC:
				double p1=0.1,p2=0.9;
				double x1 = ToolBox.percentile( intensa, p1 );
				double x2 = ToolBox.percentile( intensa, p2 );
				
				double a = (Math.log( p1/(1.0-p1) ) - Math.log( p2/(1.0-p2) ))/(x1-x2);
				double b = Math.log( p1/(1.0-p1) ) - a*x1;
				
				for(int i=0;i<weights.length;i++){
					weights[i] = 1.0/(1.0 + Math.exp( -(a*intensa[i] + b) ));
				}
				break;
			case LOG_LINEAR:
				min = Math.log(min);
				max = Math.log(max)-min;
				for(int i=0;i<intensa.length;i++){
					weights[i] = aPrioriWeight+(1-2*aPrioriWeight)*(Math.log(intensa[i]) - min)/max;
				}
				break;
			case LINEAR_jens:
				max-=min;
				for(int i=0;i<intensa.length;i++){
					weights[i] = (intensa[i] - min)/max;
				}
				break;
			case LINEAR:
				double distAbove = max - thresh;
				double distBelow = thresh - min;
				for(int i=0;i<intensa.length;i++){
					if(intensa[i] >= thresh){
						temp = intensa[i] - thresh;
						weights[i] = 0.5 + 0.5*temp/distAbove;
					}else{
						temp = thresh - intensa[i];
						weights[i] = 0.5 - 0.5*temp/distBelow;
					}
				}
				break;
			case RANK_OUTLIER_REPEATED:
			case RANK_OUTLIER:
				boolean repeat = interpol == RANK_OUTLIER_REPEATED;
				int n1, n2, idx1 = 0, idx2 = -1000;
				do {
					if( idx2 < 0 ) {
						n1 = 0;
						n2 = help.length;
					} else {
						n1 = idx1;
						n2 = idx2;
					}
					
					q25 = help[n1+(int)(0.25*(n2-n1))];
					q75 = help[n1+(int)(0.75*(n2-n1))];
					/*
					thresh = q25 - 3*(q75-q25);
					idx1 = Arrays.binarySearch( help, thresh );
					if( idx1 < 0 ) {
						idx1 = -idx1 + 1;
					}
					*/
					thresh = q75 + 1.5*(q75-q25);
					idx2 = Arrays.binarySearch( help, thresh );
					if( idx2 < 0 ) {
						idx2 = -idx2 + 1;
					}
					//System.out.println( help[n1] + "\t" + thresh + "\t" + help[n2-1] + "\t" + idx1 + "\t" + n2 );
				} while( repeat && idx2 < n2 );
				aPrioriWeight = 1 - idx2 / (double) intensa.length;
				System.out.println( aPrioriWeight );
			case RANK_LOG:
				rank = /*ToolBox.rank(intensa,true);*/ToolBox.rank( intensa, TiedRanks.SPORTS );
				for(int i=0;i<intensa.length;i++){
					if( m <= rank[i] ) {
						m = rank[i];
					}
				}
				for(int i=0;i<intensa.length;i++){
					/*
					h = 1d-(rank[i] / m);
					weights[i] = 1d/(1d + (1-h)*(1-aPrioriWeight) / (h*aPrioriWeight));//TODO changed aPrioriWeight
					weights[i] = 1d/(1d + (1-h)/h * (1-aPrioriWeight)/aPrioriWeight );
					*/
					
					h = rank[i] / m;
					weights[i] = 1d/(1d + h/(1-h) * (1-aPrioriWeight)/aPrioriWeight);
				}
				break;
			case RANK_SELEX:
				
				
				double[] weightMap = getWeightMapping( data, intensa );
				
				for(int i=0;i<weights.length;i++){
					int cyc = (int)intensa[i];
					weights[i] = weightMap[cyc];
				}
				
				//getForegroundSelexWeight( data, intensa, weights );
				
				break;
			case SIGMOID:
				b1 = (Math.log( c ) - Math.log1p( -c ))/(1.0 - max/thresh);
				a1 = b1/thresh;
				/*
				double b2 = (Math.log1p( -c ) - Math.log( c ))/(1.0 - min/thresh);
				double a2 = b2/thresh;
				if(a2 < a1){
					a1 = a2;
					b1 = b2;
				}*/
				for(int i=0;i<intensa.length;i++){
					weights[i] = 1.0/(1.0 + Math.exp( -a1*intensa[i] + b1 ));
				}
				break;
			case RANK_SIGMOID://TODO Jan
				rank = ToolBox.rank(intensa,true);
				int[] rankSorted = rank.clone();
				Arrays.sort( rankSorted );
				
				min = -rankSorted[rankSorted.length-1];
				max = -rankSorted[0];
				thresh = -rankSorted[ (int) Math.floor( (aPrioriWeight)*rankSorted.length )];
				System.out.println(min+" "+max+" "+thresh);
				
				b1 =  (Math.log1p( -c ) - Math.log( c ))/(1.0 - min/thresh);
				a1 = b1/thresh;
				/*
				double b2 = (Math.log1p( -c ) - Math.log( c ))/(1.0);
				double a2 = b2/thresh;
				if(a2 < a1){
					a1 = a2;
					b1 = b2;
				}*/
				for(int i=0;i<intensa.length;i++){
					weights[i] = 1.0/(1.0 + Math.exp( -a1*(-rank[i]) + b1 ));
				}
				break;
			default:
				/*
				AlphabetContainer contcont = new AlphabetContainer( new ContinuousAlphabet() );
				GaussianEmission em1 = new GaussianEmission( contcont );
				GaussianEmission em2 = new GaussianEmission( contcont );

				Sequence[] intenss = new Sequence[intensa.length];
				for(int i=0;i<intensa.length;i++){
					intenss[i] = new ArbitrarySequence( contcont, new double[]{intensa[i]} );
					if(intensa[i] >= thresh){
						em1.addToStatistic( true, 0, 0, 1, intenss[i] );
					}else{
						em2.addToStatistic( true, 0, 0, 1, intenss[i] );
					}
				}
				
				em1.estimateFromStatistic();
				em2.estimateFromStatistic();
				System.out.println(em1);
				System.out.println(em2);
				
				
				double logAPriori = Math.log( aPrioriWeight );
				double log1mAPriori = Math.log1p( -aPrioriWeight );
				
				double[] temp = new double[2];
				for(int i=0;i<intensa.length;i++){
					temp[0] = logAPriori + em1.getLogProbFor( true, 0, 0, intenss[i] );
					temp[1] = log1mAPriori + em2.getLogProbFor( true, 0, 0, intenss[i] );
					Normalisation.logSumNormalisation( temp );
					weights[1][0][i] = temp[0];
					weights[1][1][i] = temp[1];
				}
				*/
				throw new RuntimeException();
		}
		
		return weights;
	}
		
	
	private static double[] getForegroundSelexWeight(DataSet seqs, double[] signal, double[] weights) throws Exception {
		int max = (int) ToolBox.max( signal );
		System.out.println("max: "+max);
		double v = 0.0;
		
		while(Math.pow( v+1, v ) < max){
			v++;
		}
		System.out.println(v+": "+Math.pow( v+1, v ));
		
		
		double[] pows = new double[(int)v];
		for(int i=0;i<pows.length;i++){
			pows[i] = Math.pow( v, i );
		}
		
		LinkedList<Sequence>[] seqMap = new LinkedList[(int)v];
		for(int i=0;i<seqMap.length;i++){
			seqMap[i] = new LinkedList<Sequence>();
		}
		int[] idxs = new int[seqs.getNumberOfElements()];
		Arrays.fill( idxs, -1 );
		for(int i=0;i<seqs.getNumberOfElements();i++){
			double sig = signal[i];
			for(int j=seqMap.length-1;j>=0;j--){
				if(pows[j] <= sig){
					int idx = j;
					if(idxs[i] < 0){
						idxs[i] = idx;
					}
					seqMap[idx].add( seqs.getElementAt( i ) );
					sig -= pows[j];
				}
			}
		}
		
		double mi = Double.POSITIVE_INFINITY;
		double ma = Double.NEGATIVE_INFINITY;
		
		double[] weightMap = new double[(int)v];
		
		for(int i=0;i<seqMap.length;i++){
			if(seqMap[i].size() > 0){
				WeightedDataSetFactory fac = new WeightedDataSetFactory( SortOperation.SORT_BY_WEIGHTS, new DataSet( "", seqMap[i] ), null, 10 );
				
				double sum = fac.getSumOfWeights();
				System.out.println("num "+i+": "+fac.getNumberOfElements()+", "+fac.getSumOfWeights());
				double k = 0;
				for(int j=0;j<fac.getNumberOfElements()*0.01;j++){
					k += fac.getWeight( j );
				}
				
				
				weightMap[i] = k/sum;
				
				if(mi > weightMap[i]){
					mi = weightMap[i];
				}
				if(ma < weightMap[i]){
					ma = weightMap[i];
				}
				
			}
		}
		
		System.out.println(Arrays.toString( weightMap ));
		for(int i=0;i<weightMap.length;i++){
			if(weightMap[i] < mi){
				weightMap[i] = 0;
			}else{
				weightMap[i] = (weightMap[i] - mi)/(ma-mi)*0.9 + 0.05 ;
			}
		}
		System.out.println(Arrays.toString( weightMap ));
		
		for(int i=0;i<weights.length;i++){
			weights[i] = weightMap[idxs[i]];
		}
		
		return weights;
	}
	
	private static double[] getWeightMapping(DataSet seqs, double[] signal) throws Exception {
		int max = (int) ToolBox.max( signal );
		double[] weightMap = new double[max+1];
		
		LinkedList<Sequence>[] seqMap = new LinkedList[max+1];
		for(int i=0;i<seqMap.length;i++){
			seqMap[i] = new LinkedList<Sequence>();
		}
		
		for(int i=0;i<seqs.getNumberOfElements();i++){
			int idx = (int)signal[i];
			seqMap[idx].add( seqs.getElementAt( i ) );
		}
		
		double mi = Double.POSITIVE_INFINITY;
		double ma = Double.NEGATIVE_INFINITY;
		
		double[] ens = new double[max+1];
		double[] rels = new double[max+1];
		
		for(int i=0;i<seqMap.length;i++){
			if(seqMap[i].size() > 0){
				WeightedDataSetFactory fac = new WeightedDataSetFactory( SortOperation.SORT_BY_WEIGHTS, new DataSet( "", seqMap[i] ), null, 10 );
				
				double sum = fac.getSumOfWeights();
				System.out.println("num "+i+": "+fac.getNumberOfElements()+", "+fac.getSumOfWeights());
				double k = 0;
				for(int j=0;j<fac.getNumberOfElements()*0.01;j++){
					k += fac.getWeight( j );
				}
				
				for(int j=0;j<fac.getNumberOfElements();j++){
					ens[i] += fac.getWeight( j )/sum * Math.log( fac.getWeight( j )/sum );
				}
				
				weightMap[i] = k/sum;
				
				//double sum = fac.getSumOfWeights();
				k = 0;
				int j=0;
				while( k < sum*0.01 ){
					k += fac.getWeight( j );
					j++;
				}
			
				rels[i] = (double)j/(double)fac.getNumberOfElements();
				//weightMap[i] = (double)j/(double)fac.getNumberOfElements();
				
				if(mi > weightMap[i]){
					mi = weightMap[i];
				}
				if(ma < weightMap[i]){
					ma = weightMap[i];
				}
				
			}
		}
		System.out.println(Arrays.toString( ens ));
		System.out.println(Arrays.toString( rels ));
		System.out.println(Arrays.toString( weightMap ));
		for(int i=0;i<weightMap.length;i++){
			if(weightMap[i] < mi){
				weightMap[i] = 0;
			}else{
				weightMap[i] = (weightMap[i] - mi)/(ma-mi)*0.9 + 0.05 ;
			}
		}
		System.out.println(Arrays.toString( weightMap ));
		return weightMap;
	}
	
	
	
	public static double[] getBgWeight( double[] weights ) {
		double[] bgWeights = new double[weights.length];
		for(int i=0;i<bgWeights.length;i++){
			bgWeights[i] = 1.0 - weights[i];			
		}
		return bgWeights;
	}
}
