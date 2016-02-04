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

package de.jstacs.clustering.distances;

import de.jstacs.clustering.hierachical.ClusterTree;
import de.jstacs.data.DeBruijnGraphSequenceGenerator;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.io.ArrayHandler;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.PFMWrapperTrainSM;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.Pair;

/**
 * Helper class for comparisons of motif models based on De Bruijn sequences.
 * @author Jan Grau
 *
 */
public class DeBruijnMotifComparison {

	private static double[] getStatistics(double[] profile){
		double sum = 0.0;
		double sq = 0.0;
		for(int i=0;i<profile.length;i++){
			sum += profile[i];
			sq += profile[i]*profile[i];
		}
		return new double[]{sum,sq};
	} 
	
	private static double getCross(double[] profile1, double[] profile2, int shift){
		double cr = 0.0;
		for(int i=0;i<profile1.length;i++){
			if(i+shift>= profile2.length){
				cr += profile1[i]*profile2[i+shift-profile2.length];
			}else{
				cr += profile1[i]*profile2[i+shift];
			}
		}
		return cr;
	}
	
	
	/**
	 * Computes the correlation of the two score profiles with relative shifts of the profiles of up to <code>maxShift</code>.
	 * @param profile1 the first profile
	 * @param profile2 the second profile
	 * @param maxShift the maximum relative shift
	 * @return the maximum Pearson correlation
	 * 
	 */
	public static Pair<Integer,Double> compare(double[] profile1, double[] profile2, int maxShift){
		double[] fullStat1 = getStatistics(profile1);
		double[] fullStat2 = getStatistics(profile2);
		
		double fac = Math.sqrt( ( fullStat1[1] - 1.0/profile1.length * fullStat1[0]*fullStat1[0] ) * ( fullStat2[1] - 1.0/profile1.length * fullStat2[0]*fullStat2[0] ) );
		
		double max = Double.NEGATIVE_INFINITY;
		int maxOff = 0;
		
		for(int i=0;i<=maxShift;i++){
			
			double cr = getCross(profile1, profile2, i);
			
			double corPlus = cr - 1.0/profile1.length * fullStat1[0]*fullStat2[0];
			corPlus /= fac;
			
			cr = getCross(profile2, profile1, i);
			
			double corMinus = cr - 1.0/profile1.length * fullStat1[0]*fullStat2[0];
			corMinus /= fac;
		
			if(corPlus > max){
				max = corPlus;
				maxOff = i;
			}
			if(corMinus > max){
				max = corMinus;
				maxOff = -i;
			}
			
		}
		
		return new Pair<Integer,Double>(maxOff,max);
	}
	
	/**
	 * Returns the score profile on a De Bruin sequence for a De Bruijn sequence.
	 * @param model the model
	 * @param n the length of n-mers represented in the De Bruijn sequence
	 * @param revcom if the reverse complementary profile should be computed
	 * @param exp if the exponential scores should be returned
	 * @return the score profile
	 * @throws Exception if the score could not be computed
	 */
	public static double[][] getProfilesForMotif(StatisticalModel model, int n, boolean revcom, boolean exp) throws Exception{
		CyclicSequenceAdaptor[] ad = DeBruijnGraphSequenceGenerator.generate( (DiscreteAlphabet)model.getAlphabetContainer().getAlphabetAt( 0 ), n );
		return getProfilesForMotif( ad, model, revcom, exp );
	}
	
	/**
	 * Returns the score profile on a De Bruin sequence for a De Bruijn sequence.
	 * @param model the model
	 * @param ad the cyclic sequence(s) for which the score profile(s) should be computed
	 * @param revcom if the reverse complementary profile should be computed
	 * @param exp if the exponential scores should be returned
	 * @return the score profile
	 * @throws Exception if the score could not be computed
	 */
	public static double[][] getProfilesForMotif(CyclicSequenceAdaptor[] ad, StatisticalModel model, boolean revcom, boolean exp) throws Exception{	
		double[][] profiles = new double[ad.length][];
		
		for(int i=0;i<ad.length;i++){
			CyclicSequenceAdaptor seq = ad[i];
			if(revcom){
				seq = seq.reverseComplement();
			}
			
			int origLength = seq.getLength();
			seq = seq.getSuperSequence( seq.getLength()+model.getLength()-1 );
		
			profiles[i] = new double[origLength];
			for(int j=0;j<origLength;j++){
				if(revcom){
					if(j + model.getLength() < origLength+1){
						profiles[i][j] = model.getLogProbFor( seq, origLength-j-model.getLength()/*, origLength-j-1*/ );
					}else{
						profiles[i][j] = model.getLogProbFor( seq, seq.getLength()-j-1 /*, seq.getLength()-j-1+model.getLength()-1*/ );
					}
				}else{
					profiles[i][j] = model.getLogProbFor( seq, j/*, j+model.getLength()-1*/ );
				}
				if(exp){
					profiles[i][j] = Math.exp(profiles[i][j]);
				}
			}
		}
		return profiles;
	}
	
	/**
	 * Returns a position weight matrix (PWM) representation of the root node of the given cluster tree and
	 * also computed the relative shifts of the motifs such that they align best with the consensus motif at the root.
	 * 
	 * The integer array returned as second element of the pair is indexed in the same order as the elements
	 * of the tree elements, each element has two entries, where the first corresponds to the shift and the second is 1 for
	 * forward and -1 for reverse complementary orientation.
	 * 
	 * @param tree the tree of motifs
	 * @param n the length of n-mers used for the De Bruijn sequence for determining shifts
	 * @return the PWM and shifts (and strand orientations)
	 * @throws Exception if the score could not be computed
	 */
	public static Pair<double[][],int[][]> getClusterRepresentative(ClusterTree<StatisticalModel> tree, int n) throws Exception{
		if(tree.getNumberOfElements() == 1){
			if(tree.getClusterElements()[0] instanceof PFMWrapperTrainSM){
				return new Pair<double[][], int[][]>(((PFMWrapperTrainSM)tree.getClusterElements()[0]).getPWM(), new int[][]{{0,1}} );
			}else{
				return null;
			}
		}else{
			
			int[][] shifts = new int[tree.getNumberOfElements()][2];
			
			ClusterTree<StatisticalModel>[] subs = tree.getSubTrees();
			double[][][] reps = new double[subs.length][][];
			int[][][] prevShifts = new int[subs.length][][];
			for(int i=0;i<subs.length;i++){
				Pair<double[][],int[][]> pair = getClusterRepresentative( subs[i], n ); 
				reps[i] = pair.getFirstElement();
				prevShifts[i] = pair.getSecondElement();
			}
			double[][] rep = reps[0];
			int g = 0;
			int minPrevShift = Integer.MAX_VALUE;
			for(int i=0;i<prevShifts[0].length;i++,g++){
				shifts[g] = prevShifts[0][i].clone();
				if(shifts[g][0] < minPrevShift){
					minPrevShift = shifts[g][0];
				}
			}
			g=0;
			for(int i=0;i<prevShifts[0].length;i++,g++){
				shifts[g][0] -= minPrevShift;
			}
			
			double n1 = subs[0].getNumberOfElements();
			
			
			for(int i=1;i<reps.length;i++){
				StatisticalModel model = new PFMWrapperTrainSM( DNAAlphabetContainer.SINGLETON, null, rep, 0 );
				double[][] prof1 = getProfilesForMotif( model, n, false, false );
				StatisticalModel model2 = new PFMWrapperTrainSM( DNAAlphabetContainer.SINGLETON, null, reps[i], 0 );
				double[][] prof2 = getProfilesForMotif( model2, n, false, false );
				double[][] prof2Rc = getProfilesForMotif( model2, n, true, false );
			
				Pair<Integer,Double> fwd = compare(prof1[0], prof2[0], Math.max( rep.length - (int)Math.floor( reps[i].length ), reps[i].length - (int)Math.floor( reps.length ) ));
				
				Pair<Integer,Double> rev = compare( prof1[0], prof2Rc[0], Math.max( rep.length - (int)Math.floor( reps[i].length ), reps[i].length - (int)Math.floor( reps.length ) ) );
				int shift = fwd.getFirstElement();
				int rc = 1;
				double[][] mat = ArrayHandler.clone( reps[i] );
				if(fwd.getSecondElement() < rev.getSecondElement()){
					shift = rev.getFirstElement();
					rc = -1;
					mat = PFMComparator.getReverseComplement( DNAAlphabet.SINGLETON, mat );
				}
				int totL = shift >= 0 ? Math.max( rep.length, mat.length+shift) : Math.max( rep.length-shift, mat.length);
				double[][] com = new double[totL][rep[0].length];
				double n2 = subs[i].getNumberOfElements(); 
				for(int j=0;j<com.length;j++){
					for(int k=0;k<com[j].length;k++){
						if(shift >= 0){
							com[j][k] = ( (j < rep.length ? rep[j][k] : 0.25)*n1 + (j>= shift && (j-shift) < mat.length ? mat[j-shift][k] : 0.25)*n2 )/(n1+n2);
						}else{
							com[j][k] = ( (j >=-shift && (j+shift) < rep.length ? rep[j+shift][k] : 0.25)*n1 + (j<mat.length ? mat[j][k] : 0.25)*n2 )/(n1+n2);
						}
					}
					
				}
				n1 += n2;
				rep = com;
			
				minPrevShift = Integer.MAX_VALUE;
				for(int j=0;j<prevShifts[i].length;j++){
					prevShifts[i][j][0] *= rc;
					if(prevShifts[i][j][0] < minPrevShift){
						minPrevShift = prevShifts[i][j][0]; 
					}
					
				}
				
				
				for(int j=0;j<prevShifts[i].length;j++,g++){
					shifts[g][0] = prevShifts[i][j][0]-minPrevShift+shift;
					shifts[g][1] = prevShifts[i][j][1]*rc;
				}
			}

			return new Pair<double[][],int[][]>(rep,shifts);
		}
	}
	
}
